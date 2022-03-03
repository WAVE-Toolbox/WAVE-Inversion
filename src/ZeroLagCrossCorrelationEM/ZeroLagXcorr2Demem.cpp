#include "ZeroLagXcorr2Demem.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D emem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::ZeroLagXcorr2Demem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{
    equationType="emem"; 
    numDimension=2;
    init(ctx, dist, workflow, config);
}


template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{
    type = equationType+std::to_string(numDimension)+"D";
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrSigmaEM, ctx, dist);
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrEpsilonEM, ctx, dist);
    
    gradientDomain = config.getAndCatch("gradientDomain", 0);
    if (gradientDomain != 0) {
        IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
        dtinversion = config.get<IndexType>("DTInversion"); 
        IndexType NT = floor(ValueType(tStepEnd) / dtinversion + 0.5);
        NT = floor(ValueType(NT) / 2 + 0.5); // half recording time.
        dmemo::DistributionPtr no_dist_NT(new scai::dmemo::NoDistribution(NT));
        EXforward.setContextPtr(ctx);
        EXadjoint.setContextPtr(ctx);
        EXforward.allocate(dist, no_dist_NT);
        EXadjoint.allocate(dist, no_dist_NT);
        EYforward.setContextPtr(ctx);
        EYadjoint.setContextPtr(ctx);
        EYforward.allocate(dist, no_dist_NT);
        EYadjoint.allocate(dist, no_dist_NT);
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getContextPtr()
{
    return (xcorrEpsilonEM.getContextPtr());
}


/*! \brief override Method to write Wavefield Snapshot to file
 \param filename file name
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrSigmaEM, "xcorrSigmaEM", filename, t);
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrEpsilonEM, "xcorrEpsilonEM", filename, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrSigmaEM);
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrEpsilonEM);
}

/*! \brief Function to update the result of the zero lag cross-correlation per timestep 
 * 
\begin{equation}
\label{eqn:GradientAdjoint1}
\begin{align*}
& \nabla f_1(\varepsilon_e) = \int_0^T ( E_{1x} \pdv{E_x}{t} + E_{1y} \pdv{E_y}{t} ) dt\\
& \nabla f_1(\sigma_e) = \int_0^T ( E_{1x} E_x + E_{1y} E_y ) dt\\
\end{align*}
\end{equation}
 *
 *  Note that the forwardWavefieldDerivative is actually the derivative of the  forwardWavefield (see variable wavefieldrecordEM in IFOS.cpp).
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;    
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefield.getRefEX();
        xcorrSigmaEM += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefield.getRefEY();
        xcorrSigmaEM += temp;
        if (gradientKernel != 0) {           
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefieldDerivative.getRefEX();
        xcorrEpsilonEM += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefieldDerivative.getRefEY();
        xcorrEpsilonEM += temp;
        if (gradientKernel != 0) {     
        }
    }
}

/*! \brief Gather wavefields in the time domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::gatherWavefields(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep)
{
    scai::IndexType NT = EXforward.getNumColumns(); // half recording time.
    if (gradientKernel == 0 || decomposition == 0) {    
        if (tStep / dtinversion >= NT) {  
            EXforward.setColumn(forwardWavefield.getRefEX(), tStep / dtinversion - NT, common::BinaryOp::COPY);
            EYforward.setColumn(forwardWavefield.getRefEY(), tStep / dtinversion - NT, common::BinaryOp::COPY);
        } else {
            EXadjoint.setColumn(adjointWavefield.getRefEX(), tStep / dtinversion, common::BinaryOp::COPY);
            EYadjoint.setColumn(adjointWavefield.getRefEY(), tStep / dtinversion, common::BinaryOp::COPY);
        }
    }
}

/*! \brief Sum wavefields in the frequency domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sinFC, ValueType DT, scai::IndexType shotNumber)
{
    double start_t_shot, end_t_shot; /* For timing */
    start_t_shot = common::Walltime::get();
    
    typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;
    scai::lama::DenseMatrix<ComplexValueType> fEXforward;
    scai::lama::DenseMatrix<ComplexValueType> fEXadjoint;
    scai::lama::DenseMatrix<ComplexValueType> fEX;
    scai::lama::DenseMatrix<ComplexValueType> fEYforward;
    scai::lama::DenseMatrix<ComplexValueType> fEYadjoint;
    scai::lama::DenseMatrix<ComplexValueType> fEY;
    scai::lama::DenseVector<ComplexValueType> temp;
    scai::lama::DenseVector<ValueType> fc12Ind;
    IndexType nfc12 = 0;
    
    scai::IndexType NT = EXforward.getNumColumns(); // half recording time.
    scai::IndexType NF = sinFC.size();
    scai::IndexType nFFT = Common::calcNextPowTwo<ValueType>(NT - 1);
    ValueType df = 0.5 / (nFFT * dtinversion * DT);
    scai::lama::DenseVector<ValueType> frequencyVector = sinFC;
    ComplexValueType j(0.0, 1.0); 
    if (NF <= 1) {
        ValueType fc1 = workflow.getLowerCornerFreq();
        ValueType fc2 = workflow.getUpperCornerFreq();
        IndexType fc1Ind = ceil(fc1 / df);
        IndexType fc2Ind = ceil(fc2 / df);
        nfc12 = ceil(ValueType(fc2Ind-fc1Ind+1)/2);
        fc12Ind = lama::linearDenseVector<ValueType>(nfc12, fc1Ind, 1);
    } else {
        sinFC /= (2*df);
        sinFC += 0.5;
        fc12Ind = scai::lama::floor(sinFC);
        nfc12 = fc12Ind.size();
    }
    
    if (gradientKernel == 0 || decomposition == 0) {  
        if (gradientDomain == 1 || NF <= 1) { // FFT for all frequency samples
            fEXforward = scai::lama::cast<ComplexValueType>(EXforward);
            fEXadjoint = scai::lama::cast<ComplexValueType>(EXadjoint);
            fEXforward.resize(EXforward.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(nFFT));
            fEXadjoint.resize(EXadjoint.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(nFFT));
            fEYforward = scai::lama::cast<ComplexValueType>(EYforward);
            fEYadjoint = scai::lama::cast<ComplexValueType>(EYadjoint);
            fEYforward.resize(EYforward.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(nFFT));
            fEYadjoint.resize(EYadjoint.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(nFFT));

            scai::lama::fft<ComplexValueType>(fEXforward, 1);
            scai::lama::fft<ComplexValueType>(fEXadjoint, 1);
            scai::lama::fft<ComplexValueType>(fEYforward, 1);
            scai::lama::fft<ComplexValueType>(fEYadjoint, 1);
            scai::lama::DenseVector<ValueType> fPos = scai::lama::linearDenseVector<ValueType>(nFFT / 2 + 1, 0.0, 2*df);
            scai::lama::DenseVector<ValueType> fNeg = scai::lama::linearDenseVector<ValueType>(nFFT / 2 - 1, -(nFFT / 2 - 1) * 2*df, 2*df);
            frequencyVector.cat(fPos, fNeg);
        } else if (gradientDomain == 2 && NF > 1) { // DFT for special frequency samples.
            scai::lama::DenseMatrix<ComplexValueType> L;
            auto distNF = std::make_shared<scai::dmemo::NoDistribution>(NF);
            auto distNT = std::make_shared<scai::dmemo::NoDistribution>(NT);
            L.allocate(distNT, distNF);
            for (int tStep = 0; tStep < NT; tStep++) {
                temp = scai::lama::cast<ComplexValueType>(frequencyVector);
                temp *= -j * 2.0 * M_PI * tStep * DT;
                temp.unaryOp(temp, common::UnaryOp::EXP);
                L.setRow(temp, tStep, common::BinaryOp::COPY);
            }
            fEX = scai::lama::cast<ComplexValueType>(EXforward);
            fEXforward = fEX * L; // NR * NT * NT * NF = NR * NF
            fEX = scai::lama::cast<ComplexValueType>(EXadjoint);
            fEXadjoint = fEX * L;
            fEY = scai::lama::cast<ComplexValueType>(EYforward);
            fEYforward = fEY * L; // NR * NT * NT * NF = NR * NF
            fEY = scai::lama::cast<ComplexValueType>(EYadjoint);
            fEYadjoint = fEY * L;
        }            
        fEXadjoint.conj();
        fEYadjoint.conj();
        if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            fEX.binaryOp(fEXforward, common::BinaryOp::MULT, fEXadjoint);
            fEY.binaryOp(fEYforward, common::BinaryOp::MULT, fEYadjoint);
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fEX.getColumn(temp, fc12Ind[jf]);
                xcorrSigmaEMstep = scai::lama::real(temp);
                fEY.getColumn(temp, fc12Ind[jf]);
                xcorrEpsilonEMstep = scai::lama::real(temp);
                xcorrSigmaEMstep += xcorrEpsilonEMstep;
                if (xcorrSigmaEMstep.maxNorm() != 0)
                    xcorrSigmaEMstep *= 1.0 / xcorrSigmaEMstep.maxNorm();
                if (snapType > 0) {
                    this->writeWavefield(xcorrSigmaEMstep, "xcorrSigmaEM.step", filename, 2*fc12Ind[jf]);
                }
                xcorrSigmaEM += xcorrSigmaEMstep;
            }
        }
        if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            temp = scai::lama::cast<ComplexValueType>(frequencyVector);
            temp *= j;
            fEX.binaryOp(fEXforward, common::BinaryOp::MULT, fEXadjoint);
            fEX.scaleColumns(temp);
            fEY.binaryOp(fEYforward, common::BinaryOp::MULT, fEYadjoint);
            fEY.scaleColumns(temp);
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fEX.getColumn(temp, fc12Ind[jf]);
                xcorrEpsilonEMstep = scai::lama::real(temp);
                fEY.getColumn(temp, fc12Ind[jf]);
                xcorrSigmaEMstep = scai::lama::real(temp);
                xcorrEpsilonEMstep += xcorrSigmaEMstep;
                if (xcorrEpsilonEMstep.maxNorm() != 0)
                    xcorrEpsilonEMstep *= 1.0 / xcorrEpsilonEMstep.maxNorm();
                if (snapType > 0) {
                    this->writeWavefield(xcorrEpsilonEMstep, "xcorrEpsilonEM.step", filename, 2*fc12Ind[jf]);
                }
                xcorrEpsilonEM += xcorrEpsilonEMstep;
            }
        }
    }
    end_t_shot = common::Walltime::get();
    if (gradientDomain == 1 || NF <= 1) {
        HOST_PRINT(commShot, "Shot number " << shotNumber << ": Finish sum wavefields (FFT) in " << end_t_shot - start_t_shot << " sec.\n");
    } else if (gradientDomain == 2 && NF > 1) {
        HOST_PRINT(commShot, "Shot number " << shotNumber << ": Finish sum wavefields (DFT) in " << end_t_shot - start_t_shot << " sec.\n");
    }
}

/*! \brief Recover ZeroLagXcorr to original size
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        xcorrSigmaEM = lhs * xcorrSigmaEM;
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        xcorrEpsilonEM = lhs * xcorrEpsilonEM;
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (emem)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Getter routine for xcorrREpsilonSigmaEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getXcorrREpsilonSigmaEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrREpsilonSigmaEM in an emem modelling")
    return (xcorrREpsilonSigmaEM);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<float>;
