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
void KITGPI::SourceEstimation<ValueType>::estimateSourceSignal(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, scai::IndexType tStepEnd){
    
//     wavefield.resetWavefields();
//     
//     for (scai::IndexType tStep = 0; tStep < tStepEnd; tStep++) {
//         solver.run(receivers, sources, model, wavefield, derivatives, tStep);
//     }
// 
//     scai::lama::DenseVector<ComplexValueType> filterTmp1;
//     scai::lama::DenseVector<ComplexValueType> filterTmp2;
//     
//     addComponents(filterTmp1, receivers, receivers);
//     filterTmp1 += scai::common::Math::pow<ValueType>(waterLevel,2.0);
//     filterTmp1.unaryOp(filterTmp1, scai::common::UnaryOp::RECIPROCAL);
//     
//     addComponents(filterTmp2, receivers, receiversTrue);
//     filterTmp1 *= filterTmp2;
//     
//     filter.assignDiagonal(filterTmp1);
    applyFilter(sources);
}

/*! \brief Apply the Wiener filter to a synthetic source
 \param sources Synthetic source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::applyFilter(KITGPI::Acquisition::Sources<ValueType> &sources) {
    std::cout << sources.getSeismogramTypes().getValue(0) << std::endl;
    COMMON_THROWEXCEPTION("finished");
    
//     scai::lama::DenseMatrix<ValueType> &pSeismo = sources.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();
//     scai::lama::DenseMatrix<ComplexValueType> pSeismoTrans = scai::lama::cast<ComplexValueType>(pSeismo);
//     scai::IndexType nFFT = Common::calcNextPowTwo<ValueType>(pSeismo.getNumColumns());
//     auto colDist = std::make_shared<scai::dmemo::NoDistribution>(nFFT);
//     pSeismoTrans.resize(pSeismo.getRowDistributionPtr(),colDist);
//     scai::lama::fft<ComplexValueType>(pSeismoTrans, 0);
//     pSeismoTrans.binaryOp(pSeismoTrans, scai::common::BinaryOp::MULT, filter);
//     scai::lama::ifft<ComplexValueType>(pSeismoTrans, 0);
//     pSeismoTrans.resize(pSeismo.getRowDistributionPtr(), pSeismo.getColDistributionPtr());
//     pSeismo = scai::lama::real(pSeismoTrans);
}

/*! \brief Elementwise multiply a matrice with the conjugate complex of another matrix in frequency domain and sum the rows of the result. Fourier transform is done along the rows.
 \param prod Result
 \param A First matrix
 \param B Second matrix
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::fftMult(scai::lama::DenseVector<ComplexValueType> &prod, scai::lama::DenseMatrix<ValueType> const &A, scai::lama::DenseMatrix<ValueType> const &B) {
    scai::lama::DenseMatrix<ComplexValueType> ATmp;
    scai::lama::DenseMatrix<ComplexValueType> BTmp;
    
    ATmp = scai::lama::cast<ComplexValueType>(A);
    BTmp = scai::lama::cast<ComplexValueType>(B);
    
    scai::IndexType nFFT = Common::calcNextPowTwo<ValueType>(A.getNumColumns());
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(nFFT);
    
    ATmp.resize(A.getRowDistributionPtr(),colDist);
    BTmp.resize(A.getRowDistributionPtr(),colDist);
    
    scai::lama::fft<ComplexValueType>(ATmp, 0);
    scai::lama::fft<ComplexValueType>(BTmp, 0);
    BTmp.conj();
    ATmp.binaryOp(ATmp, scai::common::BinaryOp::MULT, BTmp);
    ATmp.reduce(prod, 1, scai::common::BinaryOp::ADD, scai::common::UnaryOp::COPY);
}

/*! \brief Multiply and sum all components of a seismogram.
 \param sum Result
 \param receiversA First receivers
 \param receiversB Second receivers
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::addComponents(scai::lama::DenseVector<ComplexValueType> sum, KITGPI::Acquisition::Receivers<ValueType> &receiversA, KITGPI::Acquisition::Receivers<ValueType> &receiversB) {
    scai::lama::DenseVector<ComplexValueType> filterTmp;
    scai::lama::DenseMatrix<ValueType> seismoA;
    scai::lama::DenseMatrix<ValueType> seismoB;

    seismoA = receiversA.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();
    seismoB = receiversB.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();
    fftMult(sum, seismoA, seismoB);
    
    seismoA = receiversA.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::VX).getData();
    seismoB = receiversB.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::VX).getData();
    fftMult(filterTmp, seismoA, seismoB);
    sum += filterTmp;
    
    seismoA = receiversA.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::VY).getData();
    seismoB = receiversB.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::VY).getData();
    fftMult(filterTmp, seismoA, seismoB);
    sum += filterTmp;
    
    seismoA = receiversA.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::VZ).getData();
    seismoB = receiversB.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::VZ).getData();
    fftMult(filterTmp, seismoA, seismoB);
    sum += filterTmp;
}

template class KITGPI::SourceEstimation<double>;
template class KITGPI::SourceEstimation<float>;