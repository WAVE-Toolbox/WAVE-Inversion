#include "SourceEstimation.hpp"

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
void KITGPI::SourceEstimation<ValueType>::estimateSourceSignal(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives){

    wavefield.resetWavefields();
    
    for (scai::IndexType tStep = 0; tStep < tStepEnd; tStep++) {
        solver.run(receivers, sources, model, wavefield, derivatives, tStep);
    }
    
    receivers.getSeismogramHandler().normalize();

    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(nFFT);
    scai::lama::DenseVector<ComplexValueType> filterTmp1(colDist, 0.0);
    scai::lama::DenseVector<ComplexValueType> filterTmp2(colDist, 0.0);
    
    addComponents(filterTmp1, receivers, receivers);
    filterTmp1 += scai::common::Math::pow<ValueType>(waterLevel,2.0);
    filterTmp1.unaryOp(filterTmp1, scai::common::UnaryOp::RECIPROCAL);
    addComponents(filterTmp2, receivers, receiversTrue);
    filterTmp1 *= filterTmp2;
    
    filter = filterTmp1; 
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
    scai::lama::DenseMatrix<ComplexValueType> seismoTransRepl;
    seismoTrans = scai::lama::cast<ComplexValueType>(seismo);
    
    // scale to power of two
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(nFFT);
    auto rowDist = std::make_shared<scai::dmemo::NoDistribution>(seismo.getNumRows());
    seismoTrans.resize(rowDist,colDist);
    
    // apply filter in frequency domain
    scai::lama::fft<ComplexValueType>(seismoTrans, 1);
    seismoTrans = scai::lama::transpose<ComplexValueType>(seismoTrans);
    seismoTrans.scaleRows(filter);
    seismoTrans = scai::lama::transpose<ComplexValueType>(seismoTrans);
    scai::lama::ifft<ComplexValueType>(seismoTrans, 1);
    
    // return to time domain
    seismoTrans.resize(seismo.getRowDistributionPtr(), seismo.getColDistributionPtr());
    seismo = scai::lama::real(seismoTrans);
    COMMON_THROWEXCEPTION("finished")
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
    BTmp.conj();
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