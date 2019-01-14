#include "SourceEstimation.hpp"

using namespace scai;

/*! \brief Initialize based on config which also initializes sourceSignalTaper
 \param config Configuration
 \param ctx Context Ptr
 \param sourceDistribution Source distribution pointer
 \param sourceSignalTaper 1D TaperPtr
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::init(Configuration::Configuration const &config, hmemo::ContextPtr ctx, dmemo::DistributionPtr sourceDistribution, std::shared_ptr<Taper::Taper<ValueType>> sourceSignalTaper)
{
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    if (config.get<bool>("useSeismogramTaper"))
        init(tStepEnd, sourceDistribution, config.get<ValueType>("waterLevel"), config.get<std::string>("seismogramTaperName"));
    else
        init(tStepEnd, sourceDistribution, config.get<ValueType>("waterLevel"));

    if (config.get<bool>("useSourceSignalTaper")) {
        sourceSignalTaper->init(std::make_shared<dmemo::NoDistribution>(tStepEnd), ctx, 1);
        if (config.get<IndexType>("sourceSignalTaperStart2") == 0 && config.get<IndexType>("sourceSignalTaperEnd2") == 0)
            sourceSignalTaper->calcCosineTaper(config.get<IndexType>("sourceSignalTaperStart1"), config.get<IndexType>("sourceSignalTaperEnd1"), 0);
        else
            sourceSignalTaper->calcCosineTaper(config.get<IndexType>("sourceSignalTaperStart1"), config.get<IndexType>("sourceSignalTaperEnd1"), config.get<IndexType>("sourceSignalTaperStart2"), config.get<IndexType>("sourceSignalTaperEnd2"), 0);
    }
}

/*! \brief Initialize based on raw variables
 \param nt Number of time steps
 \param sourceDistribution Source distribution pointer
 \param waterLvl Water level
 \param tprName if true a taper 
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::init(IndexType nt, dmemo::DistributionPtr sourceDistribution, ValueType waterLvl, std::string tprName)
{
    nFFT = Common::calcNextPowTwo<ValueType>(nt - 1);
    waterLevel = common::Math::pow<ValueType>(waterLvl, 2.0) * nFFT;
    filter.allocate(sourceDistribution, std::make_shared<dmemo::NoDistribution>(nFFT));
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
void KITGPI::SourceEstimation<ValueType>::estimateSourceSignal(KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, IndexType shotNumber)
{

    lama::DenseVector<ComplexValueType> filterTmp1(filter.getColDistributionPtr(), 0.0);
    lama::DenseVector<ComplexValueType> filterTmp2(filter.getColDistributionPtr(), 0.0);

    addComponents(filterTmp1, receivers, receivers, shotNumber);
    filterTmp1 += waterLevel * receivers.getSeismogramHandler().getNumTracesTotal();
    filterTmp1.unaryOp(filterTmp1, common::UnaryOp::RECIPROCAL);
    addComponents(filterTmp2, receivers, receiversTrue, shotNumber);
    filterTmp1 *= filterTmp2;

    filter.setRow(filterTmp1, shotNumber, common::BinaryOp::COPY);
}

/*! \brief Apply the Wiener filter to a synthetic source
 \param sources Synthetic source
 \param shotNumber Shot number of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::applyFilter(KITGPI::Acquisition::Sources<ValueType> &sources, IndexType shotNumber) const
{

    //get seismogram that corresponds to source type
    auto sourceType = Acquisition::SeismogramType(sources.getSeismogramTypes().getValue(0) - 1);
    lama::DenseMatrix<ValueType> &seismo = sources.getSeismogramHandler().getSeismogram(sourceType).getData();
    lama::DenseMatrix<ComplexValueType> seismoTrans;
    seismoTrans = lama::cast<ComplexValueType>(seismo);

    // scale to power of two
    seismoTrans.resize(seismo.getRowDistributionPtr(), filter.getColDistributionPtr());

    // apply filter in frequency domain
    lama::fft<ComplexValueType>(seismoTrans, 1);

    lama::DenseVector<ComplexValueType> filterTmp;
    filter.getRow(filterTmp, shotNumber);
    seismoTrans.scaleColumns(filterTmp);

    lama::ifft<ComplexValueType>(seismoTrans, 1);
    seismoTrans *= 1.0 / nFFT;

    // return to time domain
    seismoTrans.resize(seismo.getRowDistributionPtr(), seismo.getColDistributionPtr());
    seismo = lama::real(seismoTrans);
}

/*! \brief Correlate the rows of two matrices.
 \param prod Result
 \param A First matrix
 \param B Second matrix
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::matCorr(lama::DenseVector<ComplexValueType> &prod, lama::DenseMatrix<ValueType> const &A, lama::DenseMatrix<ValueType> const &B, IndexType iComponent)
{
    lama::DenseMatrix<ComplexValueType> ATmp;
    lama::DenseMatrix<ComplexValueType> BTmp;

    ATmp = lama::cast<ComplexValueType>(A);
    BTmp = lama::cast<ComplexValueType>(B);

    auto colDist = std::make_shared<dmemo::NoDistribution>(nFFT);
    ATmp.resize(A.getRowDistributionPtr(), colDist);
    BTmp.resize(A.getRowDistributionPtr(), colDist);

    lama::fft<ComplexValueType>(ATmp, 1);
    lama::fft<ComplexValueType>(BTmp, 1);
    ATmp.conj();
    ATmp.binaryOp(ATmp, common::BinaryOp::MULT, BTmp);

    if (useOffsetMutes)
        ATmp.scaleRows(lama::eval<lama::DenseVector<ComplexValueType>>(lama::cast<ComplexValueType>(mutes[iComponent])));

    ATmp.reduce(prod, 1, common::BinaryOp::ADD, common::UnaryOp::COPY);
}

/*! \brief Correlate and sum all components of two receivers.
 \param sum Result
 \param receiversA First receivers
 \param receiversB Second receivers
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::addComponents(lama::DenseVector<ComplexValueType> &sum, KITGPI::Acquisition::Receivers<ValueType> const &receiversA, KITGPI::Acquisition::Receivers<ValueType> const &receiversB, IndexType shotNumber)
{
    lama::DenseVector<ComplexValueType> filterTmp;

    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (receiversA.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            lama::DenseMatrix<ValueType> const &seismoA = receiversA.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            lama::DenseMatrix<ValueType> const &seismoB = receiversB.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            if (readTaper) {
                lama::DenseMatrix<ValueType> seismoA_tmp(seismoA); // tmps needed so real seismograms don't get tapered
                lama::DenseMatrix<ValueType> seismoB_tmp(seismoB);
                
                auto seismoTaper(Taper::Factory<ValueType>::Create("2D"));
                seismoTaper->init(seismoA.getRowDistributionPtr(), seismoA.getColDistributionPtr(), seismoA.getContextPtr());
                std::string taperNameTmp = taperName + ".shot_" + std::to_string(shotNumber) + "." + std::string(Acquisition::SeismogramTypeString[Acquisition::SeismogramType(iComponent)]) + ".mtx";
                seismoTaper->read(taperNameTmp, 1);
                seismoTaper->apply(seismoA_tmp);
                seismoTaper->apply(seismoB_tmp);
                
                matCorr(filterTmp, seismoA_tmp, seismoB_tmp, iComponent);
            } else
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
void KITGPI::SourceEstimation<ValueType>::calcOffsetMutes(KITGPI::Acquisition::Sources<ValueType> const &sources, KITGPI::Acquisition::Receivers<ValueType> const &receivers, ValueType maxOffset, IndexType NX, IndexType NY, IndexType NZ,KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates)
{

    lama::DenseVector<ValueType> offsets;
    lama::DenseVector<IndexType> sourceIndexVec;
    IndexType sourceIndex(0);

    bool sourceIndexSet(false);

    // get index of source position
    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (sources.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            sourceIndexVec = sources.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).get1DCoordinates();

            if (sourceIndexVec.size() > 1 || (sourceIndexSet && sourceIndexVec.getValue(0) != sourceIndex))
                COMMON_THROWEXCEPTION("offset not defined for more than one source position");

            sourceIndex = sourceIndexVec.getValue(0);
        }
    }

    // get indices of receiver position and calc mutes
    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (receivers.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            Common::calcOffsets(offsets, sourceIndex, receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).get1DCoordinates(), NX, NY, NZ,modelCoordinates);

            offsets /= -maxOffset;
            offsets.unaryOp(offsets, common::UnaryOp::CEIL);
            offsets.unaryOp(offsets, common::UnaryOp::SIGN);
            offsets += 1.0;

            mutes[iComponent] = offsets;
        }
    }
}

template class KITGPI::SourceEstimation<double>;
template class KITGPI::SourceEstimation<float>;
