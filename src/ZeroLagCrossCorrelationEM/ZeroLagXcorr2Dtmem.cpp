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
    if (workflow.getInvertForSigma() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrSigma, ctx, dist);
        if (gradientKernel != 0 && decomposition != 0) {
            this->initWavefield(xcorrSigmaSuRd, ctx, dist);
            this->initWavefield(xcorrSigmaSdRu, ctx, dist);
            this->initWavefield(xcorrSigmaSuRu, ctx, dist);
            this->initWavefield(xcorrSigmaSdRd, ctx, dist);
        }
    }
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrEpsilon, ctx, dist);
        if (gradientKernel != 0 && decomposition != 0) {
            this->initWavefield(xcorrEpsilonSuRd, ctx, dist);
            this->initWavefield(xcorrEpsilonSdRu, ctx, dist);
            this->initWavefield(xcorrEpsilonSuRu, ctx, dist);
            this->initWavefield(xcorrEpsilonSdRd, ctx, dist);
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
    return (xcorrEpsilon.getContextPtr());
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
    if (workflow.getInvertForSigma() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrSigma, "xcorrSigma", filename, t);
        if (gradientKernel == 1 && decomposition == 1) {
            this->writeWavefield(xcorrSigmaSdRu, "xcorrSigma.SdRu", filename, t);
            this->writeWavefield(xcorrSigmaSuRd, "xcorrSigma.SuRd", filename, t);
        } else if (gradientKernel == 2 && decomposition == 1) {
            this->writeWavefield(xcorrSigmaSuRu, "xcorrSigma.SuRu", filename, t);
            this->writeWavefield(xcorrSigmaSdRd, "xcorrSigma.SdRd", filename, t);
        }
    }
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
    {
        this->writeWavefield(xcorrEpsilon, "xcorrEpsilon", filename, t);
        if (gradientKernel == 1 && decomposition == 1) {
            this->writeWavefield(xcorrEpsilonSdRu, "xcorrEpsilon.SdRu", filename, t);
            this->writeWavefield(xcorrEpsilonSuRd, "xcorrEpsilon.SuRd", filename, t);
        } else if (gradientKernel == 2 && decomposition == 1) {
            this->writeWavefield(xcorrEpsilonSuRu, "xcorrEpsilon.SuRu", filename, t);
            this->writeWavefield(xcorrEpsilonSdRd, "xcorrEpsilon.SdRd", filename, t);
        }
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigma() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrSigma);
        if (gradientKernel != 0 && decomposition != 0) {
            this->resetWavefield(xcorrSigmaSuRu);
            this->resetWavefield(xcorrSigmaSdRd);
            this->resetWavefield(xcorrSigmaSuRd);
            this->resetWavefield(xcorrSigmaSdRu);
        }
    }
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrEpsilon);
        if (gradientKernel != 0 && decomposition != 0) {
            this->resetWavefield(xcorrEpsilonSuRu);
            this->resetWavefield(xcorrEpsilonSdRd);
            this->resetWavefield(xcorrEpsilonSuRd);
            this->resetWavefield(xcorrEpsilonSdRu);
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
    scai::lama::DenseVector<ValueType> temp;
    if (workflow.getInvertForSigma() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientKernel == 0 || decomposition == 0) {   
            // Born kernel or FWI kernel
            temp = adjointWavefield.getRefEZ();
            temp *= forwardWavefield.getRefEZ();
            xcorrSigma += temp;           
        } else if (gradientKernel == 1 && decomposition == 1) {    
            // migration kernel using up/down-going wavefields
            temp = adjointWavefield.getRefEZdown();
            temp *= forwardWavefield.getRefEZdown();
            xcorrSigmaSdRu += temp;           
            temp = adjointWavefield.getRefEZup();
            temp *= forwardWavefield.getRefEZup();
            xcorrSigmaSuRd += temp;      
            xcorrSigma = xcorrSigmaSdRu + xcorrSigmaSuRd;    
        } else if (gradientKernel == 2 && decomposition == 1) {    
            // tomographic kernel using up/down-going wavefields
            temp = adjointWavefield.getRefEZdown();
            temp *= forwardWavefield.getRefEZup();
            xcorrSigmaSuRu += temp;           
            temp = adjointWavefield.getRefEZup();
            temp *= forwardWavefield.getRefEZdown();
            xcorrSigmaSdRd += temp;  
            xcorrSigma = xcorrSigmaSuRu + xcorrSigmaSdRd;           
        }
    }
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientKernel == 0 || decomposition == 0) {   
            // Born kernel or FWI kernel
            temp = adjointWavefield.getRefEZ();
            temp *= forwardWavefieldDerivative.getRefEZ();
            xcorrEpsilon += temp;
        } else if (gradientKernel == 1 && decomposition == 1) {    
            // migration kernel using up/down-going wavefields   
            temp = adjointWavefield.getRefEZdown();
            temp *= forwardWavefieldDerivative.getRefEZdown();
            xcorrEpsilonSdRu += temp;           
            temp = adjointWavefield.getRefEZup();
            temp *= forwardWavefieldDerivative.getRefEZup();
            xcorrEpsilonSuRd += temp;      
            xcorrEpsilon = xcorrEpsilonSdRu + xcorrEpsilonSuRd;   
        } else if (gradientKernel == 2 && decomposition == 1) {    
            // tomographic kernel using up/down-going wavefields
            temp = adjointWavefield.getRefEZdown();
            temp *= forwardWavefieldDerivative.getRefEZup();
            xcorrEpsilonSuRu += temp;           
            temp = adjointWavefield.getRefEZup();
            temp *= forwardWavefieldDerivative.getRefEZdown();
            xcorrEpsilonSdRd += temp;  
            xcorrEpsilon = xcorrEpsilonSuRu + xcorrEpsilonSdRd;      
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sourceFC, ValueType DT, scai::IndexType shotNumber, std::vector<scai::lama::SparseVector<ValueType>> taperEncode)
{
    double start_t_shot, end_t_shot; /* For timing */
    start_t_shot = common::Walltime::get();
    
    scai::lama::DenseMatrix<ComplexValueType> ftempM;
    scai::lama::DenseVector<ComplexValueType> ftemp1;
    scai::lama::DenseVector<ComplexValueType> ftemp2;
    scai::lama::DenseVector<ValueType> xcorrSigmastep;
    scai::lama::DenseVector<ValueType> xcorrEpsilonstep;
    scai::lama::DenseVector<ValueType> fcInd; // frequency indices start from 0
    scai::lama::DenseVector<ValueType> fc12Ind; // frequency indices start from fc1
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
    fcInd = scai::lama::floor(sourceFC);
    nfc12 = fcInd.size();
    if (gradientDomain == 1) { // FFT
        fc12Ind = fcInd;
    } else if (gradientDomain == 2 || gradientDomain == 3) { // DFT or Phase sensitive detection (PSD) for special frequency samples.
        fc12Ind = scai::lama::linearDenseVector<ValueType>(nfc12, 0, 1);
    }
    
    if (gradientKernel == 0 || decomposition == 0) {  
        if (gradientDomain == 1 || NF <= 1) { // FFT for all frequency samples
            auto dist = EZforward.getRowDistributionPtr();
            auto distNFFT = std::make_shared<scai::dmemo::NoDistribution>(nFFT);
            fEZforward = scai::lama::cast<ComplexValueType>(EZforward);
            fEZadjoint = scai::lama::cast<ComplexValueType>(EZadjoint);
            fEZforward.resize(dist, distNFFT);
            fEZadjoint.resize(dist, distNFFT);
            scai::lama::fft<ComplexValueType>(fEZforward, 1);            
            scai::lama::fft<ComplexValueType>(fEZadjoint, 1);
        } else if (gradientDomain == 2 && NF > 1) { // DFT for special frequency samples.
            scai::lama::DenseMatrix<ComplexValueType> L;
            auto distNF = std::make_shared<scai::dmemo::NoDistribution>(NF);
            auto distNT = std::make_shared<scai::dmemo::NoDistribution>(NT);
            L.allocate(distNT, distNF);
            for (int tStep = 0; tStep < NT; tStep++) {
                ftemp1 = scai::lama::cast<ComplexValueType>(omega);
                ftemp1 *= -j * tStep * DT * workflow.skipDT;
                ftemp1.unaryOp(ftemp1, common::UnaryOp::EXP);
                L.setRow(ftemp1, tStep, common::BinaryOp::COPY);
            }
            ftempM = scai::lama::cast<ComplexValueType>(EZforward);
            fEZforward = ftempM * L; // NR * NT * NT * NF = NR * NF
            ftempM = scai::lama::cast<ComplexValueType>(EZadjoint);
            fEZadjoint = ftempM * L;
        }            
        
        IndexType frequencySkip = 1;
        if (snapType > 0 && useSourceEncode != 0)
            frequencySkip = 2;
        
        if (workflow.getInvertForSigma() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fEZforward.getColumn(ftemp1, fc12Ind[jf]);
                fEZadjoint.getColumn(ftemp2, fc12Ind[jf]);                
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);
                ftemp1 *= ftemp2;
                
                xcorrSigmastep = scai::lama::real(ftemp1);
                
                if (taperEncode.size() == 1) {
                    xcorrSigmastep *= taperEncode[0];
                } else if (taperEncode.size() > 1) {
                    xcorrSigmastep *= taperEncode[jf];
                }
                if (normalizeGradient && useSourceEncode != 0 && xcorrSigmastep.maxNorm() != 0)
                    xcorrSigmastep *= 1.0 / xcorrSigmastep.maxNorm();
                xcorrSigmastep *= weightingFreq[fcInd[jf]];
                if (snapType > 0 && workflow.workflowStage == 0 && workflow.iteration == 0) {
                    this->writeWavefield(xcorrSigmastep, "xcorrSigma.step", filename, frequencySkip*fcInd[jf]);
                }
                xcorrSigma += xcorrSigmastep;
            }
        }
        if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fEZforward.getColumn(ftemp1, fc12Ind[jf]);
                fEZadjoint.getColumn(ftemp2, fc12Ind[jf]);                
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);
                ftemp1 *= ftemp2;
                
                ftemp1 *= j * omega[jf];
                xcorrEpsilonstep = scai::lama::real(ftemp1);
                
                if (taperEncode.size() == 1) {
                    xcorrEpsilonstep *= taperEncode[0];
                } else if (taperEncode.size() > 1) {
                    xcorrEpsilonstep *= taperEncode[jf];
                }
                if (normalizeGradient && useSourceEncode != 0 && xcorrEpsilonstep.maxNorm() != 0)
                    xcorrEpsilonstep *= 1.0 / xcorrEpsilonstep.maxNorm();
                xcorrEpsilonstep *= weightingFreq[fcInd[jf]];
                if (snapType > 0 && useSourceEncode == 0) {
                    this->writeWavefield(xcorrEpsilonstep, "xcorrEpsilon.step", filename, frequencySkip*fcInd[jf]);
                }
                xcorrEpsilon += xcorrEpsilonstep;
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
    if (workflow.getInvertForSigma() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        xcorrSigma = lhs * xcorrSigma;
    }
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        xcorrEpsilon = lhs * xcorrEpsilon;
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

//! \brief Getter routine for xcorrREpsilonSigma
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getXcorrREpsilonSigma() const
{
    COMMON_THROWEXCEPTION("There is no xcorrREpsilonSigma in an tmem modelling")
    return (xcorrREpsilonSigma);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<float>;
