#include "SourceEstimation.hpp"

/*! \brief Initialize
 \param nt Number of time steps
 \param sourceDistribution Source distribution pointer
 \param waterLvl Water level
 \param tprName if true a taper 
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::init(scai::IndexType nt, scai::dmemo::DistributionPtr sourceDistribution, ValueType waterLvl, std::string tprName)
{
    nFFT = Common::calcNextPowTwo<ValueType>(nt - 1);
    waterLevel = scai::common::Math::pow<ValueType>(waterLvl, 2.0) * nFFT;
    filter.allocate(sourceDistribution,std::make_shared<scai::dmemo::NoDistribution>(nFFT));
    if (!tprName.empty()) {
        readTaper = true;
        taperName = tprName;
    }
}

/*! \brief Calculate the Wiener filter
 \param receivers Synthetic receivers
 \param receiversTrue Observed receivers
 \param shotNumber Shot number of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::estimateSourceSignal(KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, scai::IndexType shotNumber)
{

    scai::lama::DenseVector<ComplexValueType> filterTmp1(filter.getColDistributionPtr(), 0.0);
    scai::lama::DenseVector<ComplexValueType> filterTmp2(filter.getColDistributionPtr(), 0.0);

    addComponents(filterTmp1, receivers, receivers, shotNumber);
    filterTmp1 += waterLevel * receivers.getSeismogramHandler().getNumTracesTotal();
    filterTmp1.unaryOp(filterTmp1, scai::common::UnaryOp::RECIPROCAL);
    addComponents(filterTmp2, receivers, receiversTrue, shotNumber);
    filterTmp1 *= filterTmp2;

    filter.setRow(filterTmp1, shotNumber, scai::common::BinaryOp::COPY);
}

/*! \brief Apply the Wiener filter to a synthetic source
 \param sources Synthetic source
 \param shotNumber Shot number of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::applyFilter(KITGPI::Acquisition::Sources<ValueType> &sources, scai::IndexType shotNumber) const
{

    //get seismogram that corresponds to source type
    auto sourceType = Acquisition::SeismogramType(sources.getSeismogramTypes().getValue(0) - 1);
    scai::lama::DenseMatrix<ValueType> &seismo = sources.getSeismogramHandler().getSeismogram(sourceType).getData();
    scai::lama::DenseMatrix<ComplexValueType> seismoTrans;
    seismoTrans = scai::lama::cast<ComplexValueType>(seismo);

    // scale to power of two
    seismoTrans.resize(seismo.getRowDistributionPtr(), filter.getColDistributionPtr());

    // apply filter in frequency domain
    scai::lama::fft<ComplexValueType>(seismoTrans, 1);
    
    scai::lama::DenseVector<ComplexValueType> filterTmp;
    filter.getRow(filterTmp, shotNumber);
    seismoTrans.scaleColumns(filterTmp);
    
    scai::lama::ifft<ComplexValueType>(seismoTrans, 1);
    seismoTrans *= 1.0 / nFFT;

    // return to time domain
    seismoTrans.resize(seismo.getRowDistributionPtr(), seismo.getColDistributionPtr());
    seismo = scai::lama::real(seismoTrans);
}

/*! \brief Correlate the rows of two matrices.
 \param prod Result
 \param A First matrix
 \param B Second matrix
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::matCorr(scai::lama::DenseVector<ComplexValueType> &prod, scai::lama::DenseMatrix<ValueType> const &A, scai::lama::DenseMatrix<ValueType> const &B, scai::IndexType iComponent)
{
    scai::lama::DenseMatrix<ComplexValueType> ATmp;
    scai::lama::DenseMatrix<ComplexValueType> BTmp;

    ATmp = scai::lama::cast<ComplexValueType>(A);
    BTmp = scai::lama::cast<ComplexValueType>(B);

    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(nFFT);
    ATmp.resize(A.getRowDistributionPtr(), colDist);
    BTmp.resize(A.getRowDistributionPtr(), colDist);

    scai::lama::fft<ComplexValueType>(ATmp, 1);
    scai::lama::fft<ComplexValueType>(BTmp, 1);
    ATmp.conj();
    ATmp.binaryOp(ATmp, scai::common::BinaryOp::MULT, BTmp);
    
    if (useOffsetMutes)
        ATmp.scaleRows(scai::lama::eval<scai::lama::DenseVector<ComplexValueType>>(scai::lama::cast<ComplexValueType>(mutes[iComponent])));
    
    ATmp.reduce(prod, 1, scai::common::BinaryOp::ADD, scai::common::UnaryOp::COPY);
}

/*! \brief Correlate and sum all components of two receivers.
 \param sum Result
 \param receiversA First receivers
 \param receiversB Second receivers
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::addComponents(scai::lama::DenseVector<ComplexValueType> &sum, KITGPI::Acquisition::Receivers<ValueType> const &receiversA, KITGPI::Acquisition::Receivers<ValueType> const &receiversB, scai::IndexType shotNumber)
{
    scai::lama::DenseVector<ComplexValueType> filterTmp;

    for (scai::IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (receiversA.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            scai::lama::DenseMatrix<ValueType> const &seismoA = receiversA.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            scai::lama::DenseMatrix<ValueType> const &seismoB = receiversB.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            if (readTaper) {
                scai::lama::DenseMatrix<ValueType> seismoA_tmp(seismoA);
                scai::lama::DenseMatrix<ValueType> seismoB_tmp(seismoA);
                Taper::Taper<ValueType>::TaperPtr seismoTaper(Taper::Factory<ValueType>::Create("2D"));
                seismoTaper->init(seismoA.getRowDistributionPtr(), seismoA.getColDistributionPtr(), seismoA.getRowDistributionPtr()->getContextPtr());
                std::string taperNameTmp = taperName + ".shot_" + std::to_string(shotNumber)+ "." + std::string(Acquisition::SeismogramTypeString[Acquisition::SeismogramType(iComponent)]) + ".mtx";
                seismoTaper->read(taperName);
                seismoTaper->apply(seismoA_tmp);
                seismoTaper->apply(seismoB_tmp);
                matCorr(filterTmp, seismoA_tmp, seismoB_tmp, iComponent);
            }
            else
                matCorr(filterTmp, seismoA, seismoB, iComponent);
            sum += filterTmp;
        }
    }
}

/*! \brief Sets mutes vector needed to exclude offsets greater than a threshold from the source estimation
 \param sources Sources
 \param receivers Receivers
 \param maxOffset Offset threshold
 \param NX Number of grid points in x-direction
 \param NY Number of grid points in y-direction
 \param NZ Number of grid points in z-direction
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::calcOffsetMutes(KITGPI::Acquisition::Sources<ValueType> const &sources, KITGPI::Acquisition::Receivers<ValueType> const &receivers, ValueType maxOffset, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ) {
    
    scai::lama::DenseVector<ValueType> offsets;
    scai::lama::DenseVector<scai::IndexType> sourceIndexVec;
    scai::IndexType sourceIndex(0);
    
    bool sourceIndexSet(false);
    
    // get index of source position
    for (scai::IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (sources.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            sourceIndexVec = sources.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getCoordinates();
            
            if (sourceIndexVec.size() > 1 || (sourceIndexSet && sourceIndexVec.getValue(0) != sourceIndex))
                COMMON_THROWEXCEPTION("offset not defined for more than one source position");

            sourceIndex = sourceIndexVec.getValue(0);  
        }
    }
    
    // get indices of receiver position and calc mutes
    for (scai::IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (receivers.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            Common::calcOffsets(offsets, sourceIndex, receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getCoordinates(), NX, NY, NZ);
            
            offsets /= -maxOffset;
            offsets.unaryOp(offsets, scai::common::UnaryOp::CEIL);
            offsets.unaryOp(offsets, scai::common::UnaryOp::SIGN);
            offsets += 1.0;
            
            mutes[iComponent] = offsets;
        }
    }
}

template class KITGPI::SourceEstimation<double>;
template class KITGPI::SourceEstimation<float>;