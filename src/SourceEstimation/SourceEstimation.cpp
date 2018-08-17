#include "SourceEstimation.hpp"

/*! \brief default constructor
 */
template <typename ValueType>
KITGPI::SourceEstimation<ValueType>::SourceEstimation(scai::IndexType nt) : waterLevel(0.0), nFFT(Common::calcNextPowTwo<ValueType>(nt-1)), filter(std::make_shared<scai::dmemo::NoDistribution>(nFFT), 0.0){
}

/*! \brief Constructor which sets the water level and filter length
 \param waterLevel water level
 \param nt last time step + 1
 */
template <typename ValueType>
KITGPI::SourceEstimation<ValueType>::SourceEstimation(ValueType waterLvl,scai::IndexType nt) : nFFT(Common::calcNextPowTwo<ValueType>(nt-1)), filter(std::make_shared<scai::dmemo::NoDistribution>(nFFT), 0.0) {
    waterLevel = scai::common::Math::pow<ValueType>(waterLvl,2.0) / nFFT;
}


/*! \brief setter function for the water level
 \param waterLevel water level
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::setWaterLevel(ValueType waterLvl){
    waterLevel = scai::common::Math::pow<ValueType>(waterLvl,2.0) / nFFT;
}


/*! \brief Calculate the Wiener filter
 \param solver Forward solver
 \param receivers Synthetic receivers
 \param receiversTrue Observed receivers
 \param sources Synthetic source
 \param model Model
 \param wavefield Wavefield
 \param derivatives Derivatives
 \param tStepEnd Number of time steps
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::estimateSourceSignal(KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Acquisition::Sources<ValueType> &sources){

    scai::lama::DenseVector<ComplexValueType> filterTmp(filter.getDistributionPtr(), 0.0);
    
    addComponents(filter, receivers, receivers);
    filter += waterLevel;
    filter.unaryOp(filter, scai::common::UnaryOp::RECIPROCAL);
    addComponents(filterTmp, receivers, receiversTrue);
    filter *= filterTmp;

    applyFilter(sources);
}

/*! \brief Apply the Wiener filter to a synthetic source
 \param sources Synthetic source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::applyFilter(KITGPI::Acquisition::Sources<ValueType> &sources) {
    
    //get seismogram that corresponds to source type
    auto sourceType = Acquisition::SeismogramType(sources.getSeismogramTypes().getValue(0)-1);
    scai::lama::DenseMatrix<ValueType> &seismo = sources.getSeismogramHandler().getSeismogram(sourceType).getData();
    scai::lama::DenseMatrix<ComplexValueType> seismoTrans;
    seismoTrans = scai::lama::cast<ComplexValueType>(seismo);
    
    // scale to power of two
    seismoTrans.resize(seismo.getRowDistributionPtr(),filter.getDistributionPtr());
    
    // apply filter in frequency domain
    scai::lama::fft<ComplexValueType>(seismoTrans, 1);
    seismoTrans.scaleColumns(filter);
    scai::lama::ifft<ComplexValueType>(seismoTrans, 1);
    seismoTrans *= 1.0/nFFT;
    
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
void KITGPI::SourceEstimation<ValueType>::matCorr(scai::lama::DenseVector<ComplexValueType> &prod, scai::lama::DenseMatrix<ValueType> const &A, scai::lama::DenseMatrix<ValueType> const &B) {
    scai::lama::DenseMatrix<ComplexValueType> ATmp;
    scai::lama::DenseMatrix<ComplexValueType> BTmp;
            
    ATmp = scai::lama::cast<ComplexValueType>(A);
    BTmp = scai::lama::cast<ComplexValueType>(B);
    
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(nFFT);
    ATmp.resize(A.getRowDistributionPtr(),colDist);
    BTmp.resize(A.getRowDistributionPtr(),colDist);
        
    scai::lama::fft<ComplexValueType>(ATmp, 1);
    scai::lama::fft<ComplexValueType>(BTmp, 1);
    ATmp.conj();
    ATmp.binaryOp(ATmp, scai::common::BinaryOp::MULT, BTmp);
    ATmp.reduce(prod, 1, scai::common::BinaryOp::ADD, scai::common::UnaryOp::COPY);
}

/*! \brief Correlate and sum all components of two receivers.
 \param sum Result
 \param receiversA First receivers
 \param receiversB Second receivers
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::addComponents(scai::lama::DenseVector<ComplexValueType> &sum, KITGPI::Acquisition::Receivers<ValueType> const &receiversA, KITGPI::Acquisition::Receivers<ValueType> const &receiversB) {
    scai::lama::DenseVector<ComplexValueType> filterTmp;

    for (scai::IndexType iComponent = 0; iComponent < 4; iComponent++) {
        if (receiversA.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            scai::lama::DenseMatrix<ValueType> const &seismoA = receiversA.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            scai::lama::DenseMatrix<ValueType> const &seismoB = receiversB.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            matCorr(filterTmp, seismoA, seismoB);
            sum += filterTmp;
        }
    }
}

template class KITGPI::SourceEstimation<double>;
template class KITGPI::SourceEstimation<float>;