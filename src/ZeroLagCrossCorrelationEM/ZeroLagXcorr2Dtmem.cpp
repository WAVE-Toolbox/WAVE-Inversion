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
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::ZeroLagXcorr2Dtmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot)
{
    equationType="tmem"; 
    numDimension=2;
    init(ctx, dist, workflow, config, numShotPerSuperShot);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot)
{
    type = equationType+std::to_string(numDimension)+"D";
    gradientKernel = config.getAndCatch("gradientKernel", 0);
    decomposition = config.getAndCatch("decomposition", 0);
    useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    gradientDomain = config.getAndCatch("gradientDomain", 0);
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
    if (gradientDomain != 0) {
        IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
        NT = floor(ValueType(tStepEnd) / workflow.skipDT) + 1;
        if (useSourceEncode != 0)
            NT = floor(ValueType(NT) / 2 + 0.5); // half recording time.
        if (gradientDomain == 1 || gradientDomain == 2) {
            dmemo::DistributionPtr no_dist_NT(new scai::dmemo::NoDistribution(NT));
            EZforward.setContextPtr(ctx);
            EZadjoint.setContextPtr(ctx);
            EZforward.allocate(dist, no_dist_NT);
            EZadjoint.allocate(dist, no_dist_NT);
        } else if (gradientDomain == 3) {
            dmemo::DistributionPtr no_dist_NS(new scai::dmemo::NoDistribution(numShotPerSuperShot));
            fEZforward.setContextPtr(ctx);
            fEZadjoint.setContextPtr(ctx);
            fEZforward.allocate(dist, no_dist_NS);
            fEZadjoint.allocate(dist, no_dist_NS);
        }
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
    if (gradientDomain == 1 || gradientDomain == 2) {
        EZforward.scale(0.0);
        EZadjoint.scale(0.0);
    } else if (gradientDomain == 3) {
        fEZforward.scale(0.0);
        fEZadjoint.scale(0.0);
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::gatherWavefields(Wavefields::Wavefields<ValueType> &wavefields, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint)
{
    IndexType T0 = 0;
    if (isAdjoint || useSourceEncode != 0)
        T0 = NT; // half recording time.
    if (gradientKernel == 0 || decomposition == 0) {   
        if (gradientDomain == 1 || gradientDomain == 2) {
            if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                EZforward.setColumn(wavefields.getRefEZ(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
            } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                EZadjoint.setColumn(wavefields.getRefEZ(), tStep / workflow.skipDT, common::BinaryOp::COPY);
            }
        } else if (gradientDomain == 3) {
            scai::lama::DenseVector<ComplexValueType> temp;
            scai::lama::DenseVector<ValueType> omega(wavefields.getRefEZ().getDistributionPtr(), 2.0 * M_PI * tStep * DT);
            IndexType NF = sourceFC.size(); 
            ComplexValueType j(0.0, 1.0); 
            if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                for (IndexType jf = 0; jf < NF; jf++) {
                    temp = scai::lama::cast<ComplexValueType>(omega);
                    temp *= -j * sourceFC[jf];
                    temp.unaryOp(temp, common::UnaryOp::EXP);
                    temp *= scai::lama::cast<ComplexValueType>(wavefields.getRefEZ());
                    fEZforward.setColumn(temp, jf, common::BinaryOp::ADD);
                }
            } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                for (IndexType jf = 0; jf < NF; jf++) {
                    temp = scai::lama::cast<ComplexValueType>(omega);
                    temp *= -j * sourceFC[jf];
                    temp.unaryOp(temp, common::UnaryOp::EXP);
                    temp *= scai::lama::cast<ComplexValueType>(wavefields.getRefEZ());
                    fEZadjoint.setColumn(temp, jf, common::BinaryOp::ADD);
                }
            }
        }
    } 
}

/*! \brief Sum wavefields in the frequency domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sourceFC, ValueType DT, scai::IndexType shotNumber)
{
    double start_t_shot, end_t_shot; /* For timing */
    start_t_shot = common::Walltime::get();
    
    scai::lama::DenseMatrix<ComplexValueType> fEZ;
    scai::lama::DenseVector<ComplexValueType> temp1;
    scai::lama::DenseVector<ComplexValueType> temp2;
    scai::lama::DenseVector<ValueType> fc12Ind;
    IndexType nfc12 = 0;
    
    scai::IndexType NF = sourceFC.size();
    scai::IndexType nFFT = Common::calcNextPowTwo<ValueType>(NT - 1);
    ValueType df = 1.0 / (nFFT * workflow.skipDT * DT);
    scai::lama::DenseVector<ValueType> omega;
    scai::lama::DenseVector<ValueType> weightingFreq = workflow.getWeightingFreq();
    ComplexValueType j(0.0, 1.0); 
    
    omega = sourceFC * 2.0 * M_PI;
    sourceFC /= df;
    sourceFC += 0.5;
    fc12Ind = scai::lama::floor(sourceFC);
    nfc12 = fc12Ind.size();
    
    if (gradientKernel == 0 || decomposition == 0) {  
        if (gradientDomain == 1 || NF <= 1) { // FFT for all frequency samples
            fEZforward = scai::lama::cast<ComplexValueType>(EZforward);
            fEZadjoint = scai::lama::cast<ComplexValueType>(EZadjoint);
            fEZforward.resize(EZforward.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(nFFT));
            fEZadjoint.resize(EZadjoint.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(nFFT));

            scai::lama::fft<ComplexValueType>(fEZforward, 1);
            scai::lama::fft<ComplexValueType>(fEZadjoint, 1);
        } else if (gradientDomain == 2 && NF > 1) { // DFT for special frequency samples.
            scai::lama::DenseMatrix<ComplexValueType> L;
            auto distNF = std::make_shared<scai::dmemo::NoDistribution>(NF);
            auto distNT = std::make_shared<scai::dmemo::NoDistribution>(NT);
            L.allocate(distNT, distNF);
            for (int tStep = 0; tStep < NT; tStep++) {
                temp1 = scai::lama::cast<ComplexValueType>(omega);
                temp1 *= -j * tStep * DT * workflow.skipDT;
                temp1.unaryOp(temp1, common::UnaryOp::EXP);
                L.setRow(temp1, tStep, common::BinaryOp::COPY);
            }
            fEZ = scai::lama::cast<ComplexValueType>(EZforward);
            fEZforward = fEZ * L; // NR * NT * NT * NF = NR * NF
            fEZ = scai::lama::cast<ComplexValueType>(EZadjoint);
            fEZadjoint = fEZ * L;
        }            
        
        IndexType frequencySkip = 1;
        if (snapType > 0 && useSourceEncode != 0)
            frequencySkip = 2;
        
        if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {
                if (gradientDomain == 1) { // FFT
                    fEZforward.getColumn(temp1, fc12Ind[jf]);
                    fEZadjoint.getColumn(temp2, fc12Ind[jf]);
                } else if (gradientDomain == 2 || gradientDomain == 3) { // DFT or Phase sensitive detection (PSD) for special frequency samples.
                    fEZforward.getColumn(temp1, jf);
                    fEZadjoint.getColumn(temp2, jf);
                }
                temp2.unaryOp(temp2, common::UnaryOp::CONJ);
                temp1 *= temp2;
                xcorrSigmaEMstep = scai::lama::real(temp1);
                if (normalizeGradient && useSourceEncode != 0 && xcorrSigmaEMstep.maxNorm() != 0)
                    xcorrSigmaEMstep *= 1.0 / xcorrSigmaEMstep.maxNorm();
                xcorrSigmaEMstep *= weightingFreq[fc12Ind[jf]];
                if (snapType > 0) {
                    this->writeWavefield(xcorrSigmaEMstep, "xcorrSigmaEM.step", filename, frequencySkip*fc12Ind[jf]);
                }
                xcorrSigmaEM += xcorrSigmaEMstep;
            }
        }
        if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {
                if (gradientDomain == 1) { // FFT
                    fEZforward.getColumn(temp1, fc12Ind[jf]);
                    fEZadjoint.getColumn(temp2, fc12Ind[jf]);
                } else if (gradientDomain == 2 || gradientDomain == 3) { // DFT
                    fEZforward.getColumn(temp1, jf);
                    fEZadjoint.getColumn(temp2, jf);
                }
                temp2.unaryOp(temp2, common::UnaryOp::CONJ);
                temp1 *= temp2;
                temp1 *= j * omega[jf];
                xcorrEpsilonEMstep = scai::lama::real(temp1);
                if (normalizeGradient && useSourceEncode != 0 && xcorrEpsilonEMstep.maxNorm() != 0)
                    xcorrEpsilonEMstep *= 1.0 / xcorrEpsilonEMstep.maxNorm();
                xcorrEpsilonEMstep *= weightingFreq[fc12Ind[jf]];
                if (snapType > 0 && useSourceEncode == 0) {
                    this->writeWavefield(xcorrEpsilonEMstep, "xcorrEpsilonEM.step", filename, frequencySkip*fc12Ind[jf]);
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
    } else if (gradientDomain == 3 && NF > 1) {
        HOST_PRINT(commShot, "Shot number " << shotNumber << ": Finish sum wavefields (PSD) in " << end_t_shot - start_t_shot << " sec.\n");
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
