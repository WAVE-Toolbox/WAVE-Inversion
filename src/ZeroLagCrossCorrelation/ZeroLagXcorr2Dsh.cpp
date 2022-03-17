#include "ZeroLagXcorr2Dsh.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D sh wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::ZeroLagXcorr2Dsh(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot)
{
    equationType="sh"; 
    numDimension=2;
    init(ctx, dist, workflow, config, numShotPerSuperShot);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot)
{
    type = equationType+std::to_string(numDimension)+"D";
    gradientKernel = config.getAndCatch("gradientKernel", 0);
    decomposition = config.getAndCatch("decomposition", 0);
    useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    gradientDomain = config.getAndCatch("gradientDomain", 0);
    
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrRho, ctx, dist);
    }

    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrMuC, ctx, dist);
    }
    
    if (gradientDomain != 0) {
        IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
        NT = floor(ValueType(tStepEnd) / workflow.skipDT) + 1;
        if (useSourceEncode != 0)
            NT = floor(ValueType(NT) / 2 + 0.5); // half recording time.
        if (gradientDomain == 1 || gradientDomain == 2) {
            dmemo::DistributionPtr no_dist_NT(new scai::dmemo::NoDistribution(NT));
            if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                VZforward.setContextPtr(ctx);
                VZadjoint.setContextPtr(ctx);
                VZforward.allocate(dist, no_dist_NT);
                VZadjoint.allocate(dist, no_dist_NT);
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                Sxzforward.setContextPtr(ctx);
                Sxzadjoint.setContextPtr(ctx);
                Sxzforward.allocate(dist, no_dist_NT);
                Sxzadjoint.allocate(dist, no_dist_NT);
                Syzforward.setContextPtr(ctx);
                Syzadjoint.setContextPtr(ctx);
                Syzforward.allocate(dist, no_dist_NT);
                Syzadjoint.allocate(dist, no_dist_NT);
            }
        } else if (gradientDomain == 3) {
            dmemo::DistributionPtr no_dist_NS(new scai::dmemo::NoDistribution(numShotPerSuperShot));
            if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                fVZforward.setContextPtr(ctx);
                fVZadjoint.setContextPtr(ctx);
                fVZforward.allocate(dist, no_dist_NS);
                fVZadjoint.allocate(dist, no_dist_NS);
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                fSxzforward.setContextPtr(ctx);
                fSxzadjoint.setContextPtr(ctx);
                fSxzforward.allocate(dist, no_dist_NS);
                fSxzadjoint.allocate(dist, no_dist_NS);
                fSyzforward.setContextPtr(ctx);
                fSyzadjoint.setContextPtr(ctx);
                fSyzforward.allocate(dist, no_dist_NS);
                fSyzadjoint.allocate(dist, no_dist_NS);
            }
        }
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::getContextPtr()
{
    return (xcorrRho.getContextPtr());
}

/*! \brief override Method to write Wavefield Snapshot to file
 *
 *
 \param filename file name
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrRho, "xcorrRho", filename, t);
    }

    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrMuC, "xcorrMuC", filename, t);
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrRho);
        if (gradientDomain == 1 || gradientDomain == 2) {
            VZforward.scale(0.0);
            VZadjoint.scale(0.0);
        } else if (gradientDomain == 3) {
            fVZforward.scale(0.0);
            fVZadjoint.scale(0.0);
        }
    }

    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrMuC);
        if (gradientDomain == 1 || gradientDomain == 2) {
            Sxzforward.scale(0.0);
            Sxzadjoint.scale(0.0);
            Syzforward.scale(0.0);
            Syzadjoint.scale(0.0);
        } else if (gradientDomain == 3) {
            fSxzforward.scale(0.0);
            fSxzadjoint.scale(0.0);
            fSyzforward.scale(0.0);
            fSyzadjoint.scale(0.0);
        }
    }
}

/*! \brief Function to update the result of the zero lag cross-correlation per timestep 
 * 
 * The zero lag cross-correlation, \f$ X \f$, is updated with the following equations where index "forw" refers to the forward propagated wavefield and "adj" to the adjoint wavefield:
 \f{eqnarray*}
   X_{\mu} &+=& P_{\mathrm{forw}} \cdot P_{\mathrm{adj}}  \\ 
   X_{\mu} &+=& V_{x,\mathrm{forw}} \cdot V_{x,\mathrm{adj}} + V_{y,\mathrm{forw}} \cdot V_{y,\mathrm{adj}}
 \f}
 *
 * 
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefieldDerivative.getRefVZ();
        temp *= adjointWavefield.getRefVZ();
        xcorrRho += temp;
    }
    
    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefieldDerivative.getRefSxz();
        temp *= adjointWavefield.getRefSxz();
        xcorrMuC += temp;
        temp = forwardWavefieldDerivative.getRefSyz();
        temp *= adjointWavefield.getRefSyz();
        xcorrMuC += temp;
    }
}

/*! \brief Gather wavefields in the time domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::gatherWavefields(Wavefields::Wavefields<ValueType> &wavefields, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint)
{
    IndexType T0 = 0;
    if (isAdjoint || useSourceEncode != 0)
        T0 = NT; // half recording time.
    if (gradientKernel == 0 || decomposition == 0) {   
        if (gradientDomain == 1 || gradientDomain == 2) {
            if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                    VZforward.setColumn(wavefields.getRefVZ(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
                } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                    VZadjoint.setColumn(wavefields.getRefVZ(), tStep / workflow.skipDT, common::BinaryOp::COPY);
                }
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                    Sxzforward.setColumn(wavefields.getRefSxz(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
                    Syzforward.setColumn(wavefields.getRefSyz(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
                } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                    Sxzadjoint.setColumn(wavefields.getRefSxz(), tStep / workflow.skipDT, common::BinaryOp::COPY);
                    Syzadjoint.setColumn(wavefields.getRefSyz(), tStep / workflow.skipDT, common::BinaryOp::COPY);
                }
            }
        } else if (gradientDomain == 3) {
            scai::lama::DenseVector<ComplexValueType> temp;
            scai::lama::DenseVector<ComplexValueType> temp1;
            scai::lama::DenseVector<ValueType> omega(wavefields.getRefVZ().getDistributionPtr(), 2.0 * M_PI * tStep * DT);
            IndexType NF = sourceFC.size(); 
            ComplexValueType j(0.0, 1.0); 
            if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                    for (IndexType jf = 0; jf < NF; jf++) {
                        temp1 = scai::lama::cast<ComplexValueType>(omega);
                        temp1 *= -j * sourceFC[jf];
                        temp1.unaryOp(temp1, common::UnaryOp::EXP);
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefVZ());
                        fVZforward.setColumn(temp1, jf, common::BinaryOp::ADD);
                    }
                } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                    for (IndexType jf = 0; jf < NF; jf++) {
                        temp1 = scai::lama::cast<ComplexValueType>(omega);
                        temp1 *= -j * sourceFC[jf];
                        temp1.unaryOp(temp1, common::UnaryOp::EXP);
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefVZ());
                        fVZadjoint.setColumn(temp1, jf, common::BinaryOp::ADD);
                    }
                }
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                    for (IndexType jf = 0; jf < NF; jf++) {
                        temp = scai::lama::cast<ComplexValueType>(omega);
                        temp *= -j * sourceFC[jf];
                        temp.unaryOp(temp, common::UnaryOp::EXP);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefSxz());
                        fSxzforward.setColumn(temp1, jf, common::BinaryOp::ADD);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefSyz());
                        fSyzforward.setColumn(temp1, jf, common::BinaryOp::ADD);
                    }
                } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                    for (IndexType jf = 0; jf < NF; jf++) {
                        temp = scai::lama::cast<ComplexValueType>(omega);
                        temp *= -j * sourceFC[jf];
                        temp.unaryOp(temp, common::UnaryOp::EXP);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefSxz());
                        fSxzadjoint.setColumn(temp1, jf, common::BinaryOp::ADD);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefSyz());
                        fSyzadjoint.setColumn(temp1, jf, common::BinaryOp::ADD);
                    }
                }
            }
        }
    } 
}


/*! \brief Sum wavefields in the frequency domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sourceFC, ValueType DT, scai::IndexType shotNumber, std::vector<scai::lama::SparseVector<ValueType>> taperEncode)
{
    double start_t_shot, end_t_shot; /* For timing */
    start_t_shot = common::Walltime::get();
    
    scai::lama::DenseMatrix<ComplexValueType> ftempM;
    scai::lama::DenseVector<ComplexValueType> ftemp;
    scai::lama::DenseVector<ComplexValueType> ftemp1;
    scai::lama::DenseVector<ComplexValueType> ftemp2;
    scai::lama::DenseVector<ValueType> xcorrRhostep;
    scai::lama::DenseVector<ValueType> xcorrMuCstep;
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
            auto dist = VZforward.getRowDistributionPtr();
            auto distNFFT = std::make_shared<scai::dmemo::NoDistribution>(nFFT);
            if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                fVZforward = scai::lama::cast<ComplexValueType>(VZforward);
                fVZadjoint = scai::lama::cast<ComplexValueType>(VZadjoint);
                fVZforward.resize(dist, distNFFT);
                fVZadjoint.resize(dist, distNFFT);
                scai::lama::fft<ComplexValueType>(fVZforward, 1);
                scai::lama::fft<ComplexValueType>(fVZadjoint, 1);
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                fSxzforward = scai::lama::cast<ComplexValueType>(Sxzforward);
                fSxzadjoint = scai::lama::cast<ComplexValueType>(Sxzadjoint);
                fSxzforward.resize(dist, distNFFT);
                fSxzadjoint.resize(dist, distNFFT);
                scai::lama::fft<ComplexValueType>(fSxzforward, 1);
                scai::lama::fft<ComplexValueType>(fSxzadjoint, 1);

                fSyzforward = scai::lama::cast<ComplexValueType>(Syzforward);
                fSyzadjoint = scai::lama::cast<ComplexValueType>(Syzadjoint);
                fSyzforward.resize(dist, distNFFT);
                fSyzadjoint.resize(dist, distNFFT);
                scai::lama::fft<ComplexValueType>(fSyzforward, 1);
                scai::lama::fft<ComplexValueType>(fSyzadjoint, 1);
            }
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
            if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                ftempM = scai::lama::cast<ComplexValueType>(VZforward);
                fVZforward = ftempM * L; // NR * NT * NT * NF = NR * NF
                ftempM = scai::lama::cast<ComplexValueType>(VZadjoint);
                fVZadjoint = ftempM * L;
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                ftempM = scai::lama::cast<ComplexValueType>(Sxzforward);
                fSxzforward = ftempM * L; // NR * NT * NT * NF = NR * NF
                ftempM = scai::lama::cast<ComplexValueType>(Sxzadjoint);
                fSxzadjoint = ftempM * L;
                ftempM = scai::lama::cast<ComplexValueType>(Syzforward);
                fSyzforward = ftempM * L; // NR * NT * NT * NF = NR * NF
                ftempM = scai::lama::cast<ComplexValueType>(Syzadjoint);
                fSyzadjoint = ftempM * L;
            }
        }            
        
        IndexType frequencySkip = 1;
        if (snapType > 0 && useSourceEncode != 0)
            frequencySkip = 2;
        
        if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {                
                fVZforward.getColumn(ftemp1, fc12Ind[jf]);
                fVZadjoint.getColumn(ftemp2, fc12Ind[jf]);              
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);  
                ftemp1 *= ftemp2;
                
                xcorrRhostep = scai::lama::real(ftemp1);
                
                if (normalizeGradient && useSourceEncode != 0 && xcorrRhostep.maxNorm() != 0)
                    xcorrRhostep *= 1.0 / xcorrRhostep.maxNorm();
                xcorrRhostep *= weightingFreq[fcInd[jf]];
                if (snapType > 0 && workflow.workflowStage == 0 && workflow.iteration == 0) {
                    this->writeWavefield(xcorrRhostep, "xcorrRho.step", filename, frequencySkip*fcInd[jf]);
                }
                xcorrRho += xcorrRhostep;
            }
        }
        if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fSxzforward.getColumn(ftemp1, fc12Ind[jf]);
                fSxzadjoint.getColumn(ftemp2, fc12Ind[jf]);               
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);                
                ftemp = ftemp1;
                ftemp *= ftemp2;
                
                fSyzforward.getColumn(ftemp1, fc12Ind[jf]);
                fSyzadjoint.getColumn(ftemp2, fc12Ind[jf]); 
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);
                ftemp1 *= ftemp2;
                ftemp1 += ftemp;
                
                xcorrMuCstep = scai::lama::real(ftemp1);
                
                if (normalizeGradient && useSourceEncode != 0 && xcorrMuCstep.maxNorm() != 0) {
                    xcorrMuCstep *= 1.0 / xcorrMuCstep.maxNorm();
                }
                xcorrMuCstep *= weightingFreq[fcInd[jf]];
                if (snapType > 0 && workflow.workflowStage == 0 && workflow.iteration == 0) {
                    this->writeWavefield(xcorrMuCstep, "xcorrMuC.step", filename, frequencySkip*fcInd[jf]);
                }
                xcorrMuC += xcorrMuCstep;
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        xcorrMuC = lhs * xcorrMuC;
    }
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {       
        xcorrRho = lhs * xcorrRho;
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (sh)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::getEquationType() const
{
    return (equationType);
}


//! \brief Not valid in the sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::getXcorrMuA() const
{
    COMMON_THROWEXCEPTION("There is no MuA Gradient in the sh case.");
    return (xcorrMuA);
}

//! \brief Not valid in the sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::getXcorrMuB() const
{
    COMMON_THROWEXCEPTION("There is no MuB Gradient in the sh case.");
    return (xcorrMuB);
}

//! \brief Not valid in the sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::getXcorrLambda() const
{
    COMMON_THROWEXCEPTION("There is no Lambda Gradient in the sh case.");
    return (xcorrLambda);
}
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<float>;
