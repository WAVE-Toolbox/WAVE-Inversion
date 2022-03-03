#include "ZeroLagXcorr2Dtmem.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D tmem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::ZeroLagXcorr2Dtmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{
    equationType="tmem"; 
    numDimension=2;
    init(ctx, dist, workflow, config);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{
    type = equationType+std::to_string(numDimension)+"D";
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrSigmaEM, ctx, dist);
        if (gradientKernel != 0 && decomposition != 0) {
            this->initWavefield(xcorrSigmaEMSuRd, ctx, dist);
            this->initWavefield(xcorrSigmaEMSdRu, ctx, dist);
            this->initWavefield(xcorrSigmaEMSuRu, ctx, dist);
            this->initWavefield(xcorrSigmaEMSdRd, ctx, dist);
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrEpsilonEM, ctx, dist);
        if (gradientKernel != 0 && decomposition != 0) {
            this->initWavefield(xcorrEpsilonEMSuRd, ctx, dist);
            this->initWavefield(xcorrEpsilonEMSdRu, ctx, dist);
            this->initWavefield(xcorrEpsilonEMSuRu, ctx, dist);
            this->initWavefield(xcorrEpsilonEMSdRd, ctx, dist);
        }
    }
    gradientDomain = config.getAndCatch("gradientDomain", 0);
    if (gradientDomain != 0) {
        IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
        dtinversion = config.get<IndexType>("DTInversion"); 
        IndexType NT = floor(ValueType(tStepEnd) / dtinversion + 0.5);
        NT = floor(ValueType(NT) / 2 + 0.5); // half recording time.
        dmemo::DistributionPtr no_dist_NT(new scai::dmemo::NoDistribution(NT));
        EZforward.setContextPtr(ctx);
        EZadjoint.setContextPtr(ctx);
        EZforward.allocate(dist, no_dist_NT);
        EZadjoint.allocate(dist, no_dist_NT);
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getContextPtr()
{
    return (xcorrEpsilonEM.getContextPtr());
}

/*! \brief override Method to write Wavefield Snapshot to file
 *
 *
 \param filename file name
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrSigmaEM, "xcorrSigmaEM", filename, t);
        if (decomposition == 0) {
            this->writeWavefield(xcorrSigmaEMstep, "xcorrSigmaEM.step", filename, t);
        } else if (gradientKernel == 1 && decomposition == 1) {
            this->writeWavefield(xcorrSigmaEMSdRu, "xcorrSigmaEM.SdRu", filename, t);
            this->writeWavefield(xcorrSigmaEMSuRd, "xcorrSigmaEM.SuRd", filename, t);
        } else if (gradientKernel == 2 && decomposition == 1) {
            this->writeWavefield(xcorrSigmaEMSuRu, "xcorrSigmaEM.SuRu", filename, t);
            this->writeWavefield(xcorrSigmaEMSdRd, "xcorrSigmaEM.SdRd", filename, t);
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
    {
        this->writeWavefield(xcorrEpsilonEM, "xcorrEpsilonEM", filename, t);
        if (decomposition == 0) {
            this->writeWavefield(xcorrEpsilonEMstep, "xcorrEpsilonEM.step", filename, t);
        } else if (gradientKernel == 1 && decomposition == 1) {
            this->writeWavefield(xcorrEpsilonEMSdRu, "xcorrEpsilonEM.SdRu", filename, t);
            this->writeWavefield(xcorrEpsilonEMSuRd, "xcorrEpsilonEM.SuRd", filename, t);
        } else if (gradientKernel == 2 && decomposition == 1) {
            this->writeWavefield(xcorrEpsilonEMSuRu, "xcorrEpsilonEM.SuRu", filename, t);
            this->writeWavefield(xcorrEpsilonEMSdRd, "xcorrEpsilonEM.SdRd", filename, t);
        }
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrSigmaEM);
        if (gradientKernel != 0 && decomposition != 0) {
            this->resetWavefield(xcorrSigmaEMSuRu);
            this->resetWavefield(xcorrSigmaEMSdRd);
            this->resetWavefield(xcorrSigmaEMSuRd);
            this->resetWavefield(xcorrSigmaEMSdRu);
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrEpsilonEM);
        if (gradientKernel != 0 && decomposition != 0) {
            this->resetWavefield(xcorrEpsilonEMSuRu);
            this->resetWavefield(xcorrEpsilonEMSdRd);
            this->resetWavefield(xcorrEpsilonEMSuRd);
            this->resetWavefield(xcorrEpsilonEMSdRu);
        }
    }
    if (gradientDomain != 0) {
        EZforward.scale(0.0);
        EZadjoint.scale(0.0);
    }
}

/*! \brief Function to update the result of the zero lag cross-correlation per timestep 
 * 
\begin{equation}
\label{eqn:GradientAdjoint1}
\begin{align*}
& \nabla f_1(\varepsilon_e) = \int_0^T ( E_{1z} \pdv{E_z}{t} ) dt\\
& \nabla f_1(\sigma_e) = \int_0^T ( E_{1z} E_z  ) dt\\
\end{align*}
\end{equation}
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientKernel == 0 || decomposition == 0) {   
            // Born kernel or FWI kernel
            xcorrSigmaEMstep = adjointWavefield.getRefEZ();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZ();
            xcorrSigmaEM += xcorrSigmaEMstep;           
        } else if (gradientKernel == 1 && decomposition == 1) {    
            // migration kernel using up/down-going wavefields
            xcorrSigmaEMstep = adjointWavefield.getRefEZdown();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZdown();
            xcorrSigmaEMSdRu += xcorrSigmaEMstep;           
            xcorrSigmaEMstep = adjointWavefield.getRefEZup();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZup();
            xcorrSigmaEMSuRd += xcorrSigmaEMstep;      
            xcorrSigmaEM = xcorrSigmaEMSdRu + xcorrSigmaEMSuRd;    
        } else if (gradientKernel == 2 && decomposition == 1) {    
            // tomographic kernel using up/down-going wavefields
            xcorrSigmaEMstep = adjointWavefield.getRefEZdown();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZup();
            xcorrSigmaEMSuRu += xcorrSigmaEMstep;           
            xcorrSigmaEMstep = adjointWavefield.getRefEZup();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZdown();
            xcorrSigmaEMSdRd += xcorrSigmaEMstep;  
            xcorrSigmaEM = xcorrSigmaEMSuRu + xcorrSigmaEMSdRd;           
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientKernel == 0 || decomposition == 0) {   
            // Born kernel or FWI kernel
            xcorrEpsilonEMstep = adjointWavefield.getRefEZ();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZ();
            xcorrEpsilonEM += xcorrEpsilonEMstep;
        } else if (gradientKernel == 1 && decomposition == 1) {    
            // migration kernel using up/down-going wavefields   
            xcorrEpsilonEMstep = adjointWavefield.getRefEZdown();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZdown();
            xcorrEpsilonEMSdRu += xcorrEpsilonEMstep;           
            xcorrEpsilonEMstep = adjointWavefield.getRefEZup();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZup();
            xcorrEpsilonEMSuRd += xcorrEpsilonEMstep;      
            xcorrEpsilonEM = xcorrEpsilonEMSdRu + xcorrEpsilonEMSuRd;   
        } else if (gradientKernel == 2 && decomposition == 1) {    
            // tomographic kernel using up/down-going wavefields
            xcorrEpsilonEMstep = adjointWavefield.getRefEZdown();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZup();
            xcorrEpsilonEMSuRu += xcorrEpsilonEMstep;           
            xcorrEpsilonEMstep = adjointWavefield.getRefEZup();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZdown();
            xcorrEpsilonEMSdRd += xcorrEpsilonEMstep;  
            xcorrEpsilonEM = xcorrEpsilonEMSuRu + xcorrEpsilonEMSdRd;      
        }
    }
}

/*! \brief Gather wavefields in the time domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::gatherWavefields(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep)
{
    scai::IndexType NT = EZforward.getNumColumns(); // half recording time.
    if (gradientKernel == 0 || decomposition == 0) {   
        if (tStep / dtinversion >= NT) {
            EZforward.setColumn(forwardWavefield.getRefEZ(), tStep / dtinversion - NT, common::BinaryOp::COPY);
        } else {
            EZadjoint.setColumn(adjointWavefield.getRefEZ(), tStep / dtinversion, common::BinaryOp::COPY);
        }
    } 
}

/*! \brief Sum wavefields in the frequency domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sinFC, ValueType DT, scai::IndexType shotNumber)
{
    double start_t_shot, end_t_shot; /* For timing */
    start_t_shot = common::Walltime::get();
    
    typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;
    scai::lama::DenseMatrix<ComplexValueType> fEZforward;
    scai::lama::DenseMatrix<ComplexValueType> fEZadjoint;
    scai::lama::DenseMatrix<ComplexValueType> fEZ;
    scai::lama::DenseVector<ComplexValueType> temp;
    scai::lama::DenseVector<ValueType> fc12Ind;
    IndexType nfc12 = 0;
    
    scai::IndexType NT = EZforward.getNumColumns(); // half recording time.
    scai::IndexType NF = sinFC.size();
    scai::IndexType nFFT = Common::calcNextPowTwo<ValueType>(NT - 1);
    ValueType df = 0.5 / (nFFT * dtinversion * DT);
    scai::lama::DenseVector<ValueType> frequencyVector = sinFC;
    scai::lama::DenseVector<ValueType> weightingFreq = workflow.getWeightingFreq();
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
            fEZforward = scai::lama::cast<ComplexValueType>(EZforward);
            fEZadjoint = scai::lama::cast<ComplexValueType>(EZadjoint);
            fEZforward.resize(EZforward.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(nFFT));
            fEZadjoint.resize(EZadjoint.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(nFFT));

            scai::lama::fft<ComplexValueType>(fEZforward, 1);
            scai::lama::fft<ComplexValueType>(fEZadjoint, 1);
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
            fEZ = scai::lama::cast<ComplexValueType>(EZforward);
            fEZforward = fEZ * L; // NR * NT * NT * NF = NR * NF
            fEZ = scai::lama::cast<ComplexValueType>(EZadjoint);
            fEZadjoint = fEZ * L;
        }            
        fEZadjoint.conj();
        if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            fEZ.binaryOp(fEZforward, common::BinaryOp::MULT, fEZadjoint);
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fEZ.getColumn(temp, fc12Ind[jf]);
                xcorrSigmaEMstep = scai::lama::real(temp);
                if (xcorrSigmaEMstep.maxNorm() != 0)
                    xcorrSigmaEMstep *= 1.0 / xcorrSigmaEMstep.maxNorm();
                xcorrSigmaEMstep *= weightingFreq[fc12Ind[jf]];
                if (snapType > 0) {
                    this->writeWavefield(xcorrSigmaEMstep, "xcorrSigmaEM.step", filename, 2*fc12Ind[jf]);
                }
                xcorrSigmaEM += xcorrSigmaEMstep;
            }
        }
        if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            temp = scai::lama::cast<ComplexValueType>(frequencyVector);
            temp *= j;
            fEZ.binaryOp(fEZforward, common::BinaryOp::MULT, fEZadjoint);
            fEZ.scaleColumns(temp);
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fEZ.getColumn(temp, fc12Ind[jf]);
                xcorrEpsilonEMstep = scai::lama::real(temp);
                if (xcorrEpsilonEMstep.maxNorm() != 0)
                    xcorrEpsilonEMstep *= 1.0 / xcorrEpsilonEMstep.maxNorm();
                xcorrEpsilonEMstep *= weightingFreq[fc12Ind[jf]];
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Workflow::Workflow<ValueType> const &workflow)
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
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (tmem)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Getter routine for xcorrREpsilonSigmaEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getXcorrREpsilonSigmaEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrREpsilonSigmaEM in an tmem modelling")
    return (xcorrREpsilonSigmaEM);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<float>;
