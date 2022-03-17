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
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::ZeroLagXcorr2Demem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot)
{
    equationType="emem"; 
    numDimension=2;
    init(ctx, dist, workflow, config, numShotPerSuperShot);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot)
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
            EXforward.setContextPtr(ctx);
            EXadjoint.setContextPtr(ctx);
            EXforward.allocate(dist, no_dist_NT);
            EXadjoint.allocate(dist, no_dist_NT);
            EYforward.setContextPtr(ctx);
            EYadjoint.setContextPtr(ctx);
            EYforward.allocate(dist, no_dist_NT);
            EYadjoint.allocate(dist, no_dist_NT);
        } else if (gradientDomain == 3) {
            dmemo::DistributionPtr no_dist_NS(new scai::dmemo::NoDistribution(numShotPerSuperShot));
            fEXforward.setContextPtr(ctx);
            fEXadjoint.setContextPtr(ctx);
            fEXforward.allocate(dist, no_dist_NS);
            fEXadjoint.allocate(dist, no_dist_NS);
            fEYforward.setContextPtr(ctx);
            fEYadjoint.setContextPtr(ctx);
            fEYforward.allocate(dist, no_dist_NS);
            fEYadjoint.allocate(dist, no_dist_NS);
        }
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
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrSigmaEM, "xcorrSigmaEM", filename, t);
        if (gradientKernel == 1 && decomposition == 1) {
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
        if (gradientKernel == 1 && decomposition == 1) {
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
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
        EXforward.scale(0.0);
        EXadjoint.scale(0.0);
        EYforward.scale(0.0);
        EYadjoint.scale(0.0);
    } else if (gradientDomain == 3) {
        fEXforward.scale(0.0);
        fEXadjoint.scale(0.0);
        fEYforward.scale(0.0);
        fEYadjoint.scale(0.0);
    }
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
    scai::lama::DenseVector<ValueType> temp;
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientKernel == 0 || decomposition == 0) {   
            // Born kernel or FWI kernel
            temp = adjointWavefield.getRefEX();
            temp *= forwardWavefield.getRefEX();
            xcorrSigmaEM += temp;    
            temp = adjointWavefield.getRefEY();
            temp *= forwardWavefield.getRefEY();
            xcorrSigmaEM += temp;         
        } else if (gradientKernel == 1 && decomposition == 1) {    
            // migration kernel using up/down-going wavefields
            temp = adjointWavefield.getRefEXdown();
            temp *= forwardWavefield.getRefEXdown();
            xcorrSigmaEMSdRu += temp;           
            temp = adjointWavefield.getRefEYdown();
            temp *= forwardWavefield.getRefEYdown();
            xcorrSigmaEMSdRu += temp;           
            temp = adjointWavefield.getRefEXup();
            temp *= forwardWavefield.getRefEXup();
            xcorrSigmaEMSuRd += temp;      
            temp = adjointWavefield.getRefEYup();
            temp *= forwardWavefield.getRefEYup();
            xcorrSigmaEMSuRd += temp;      
            xcorrSigmaEM = xcorrSigmaEMSdRu + xcorrSigmaEMSuRd;    
        } else if (gradientKernel == 2 && decomposition == 1) {    
            // tomographic kernel using up/down-going wavefields
            temp = adjointWavefield.getRefEXdown();
            temp *= forwardWavefield.getRefEXup();
            xcorrSigmaEMSuRu += temp;           
            temp = adjointWavefield.getRefEYdown();
            temp *= forwardWavefield.getRefEYup();
            xcorrSigmaEMSuRu += temp;     
            temp = adjointWavefield.getRefEXup();
            temp *= forwardWavefield.getRefEXdown();
            xcorrSigmaEMSdRd += temp;  
            temp = adjointWavefield.getRefEYup();
            temp *= forwardWavefield.getRefEYdown();
            xcorrSigmaEMSdRd += temp;  
            xcorrSigmaEM = xcorrSigmaEMSuRu + xcorrSigmaEMSdRd;           
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientKernel == 0 || decomposition == 0) {   
            // Born kernel or FWI kernel
            temp = adjointWavefield.getRefEX();
            temp *= forwardWavefieldDerivative.getRefEX();
            xcorrEpsilonEM += temp;
            temp = adjointWavefield.getRefEY();
            temp *= forwardWavefieldDerivative.getRefEY();
            xcorrEpsilonEM += temp;
        } else if (gradientKernel == 1 && decomposition == 1) {    
            // migration kernel using up/down-going wavefields   
            temp = adjointWavefield.getRefEXdown();
            temp *= forwardWavefieldDerivative.getRefEXdown();
            xcorrEpsilonEMSdRu += temp;        
            temp = adjointWavefield.getRefEYdown();
            temp *= forwardWavefieldDerivative.getRefEYdown();
            xcorrEpsilonEMSdRu += temp;            
            temp = adjointWavefield.getRefEXup();
            temp *= forwardWavefieldDerivative.getRefEXup();
            xcorrEpsilonEMSuRd += temp;        
            temp = adjointWavefield.getRefEYup();
            temp *= forwardWavefieldDerivative.getRefEYup();
            xcorrEpsilonEMSuRd += temp;      
            xcorrEpsilonEM = xcorrEpsilonEMSdRu + xcorrEpsilonEMSuRd;   
        } else if (gradientKernel == 2 && decomposition == 1) {    
            // tomographic kernel using up/down-going wavefields
            temp = adjointWavefield.getRefEXdown();
            temp *= forwardWavefieldDerivative.getRefEXup();
            xcorrEpsilonEMSuRu += temp;   
            temp = adjointWavefield.getRefEYdown();
            temp *= forwardWavefieldDerivative.getRefEYup();
            xcorrEpsilonEMSuRu += temp;           
            temp = adjointWavefield.getRefEXup();
            temp *= forwardWavefieldDerivative.getRefEXdown();
            xcorrEpsilonEMSdRd += temp;          
            temp = adjointWavefield.getRefEYup();
            temp *= forwardWavefieldDerivative.getRefEYdown();
            xcorrEpsilonEMSdRd += temp;  
            xcorrEpsilonEM = xcorrEpsilonEMSuRu + xcorrEpsilonEMSdRd;      
        }
    }
}

/*! \brief Gather wavefields in the time domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::gatherWavefields(Wavefields::Wavefields<ValueType> &wavefields, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint)
{
    IndexType T0 = 0;
    if (isAdjoint || useSourceEncode != 0)
        T0 = NT; // half recording time.
    if (gradientKernel == 0 || decomposition == 0) {   
        if (gradientDomain == 1 || gradientDomain == 2) {
            if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                EXforward.setColumn(wavefields.getRefEX(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
                EYforward.setColumn(wavefields.getRefEY(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
            } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                EXadjoint.setColumn(wavefields.getRefEX(), tStep / workflow.skipDT, common::BinaryOp::COPY);
                EYadjoint.setColumn(wavefields.getRefEY(), tStep / workflow.skipDT, common::BinaryOp::COPY);
            }
        } else if (gradientDomain == 3) {
            scai::lama::DenseVector<ComplexValueType> temp;
            scai::lama::DenseVector<ComplexValueType> ftemp1;
            scai::lama::DenseVector<ValueType> omega(wavefields.getRefEX().getDistributionPtr(), 2.0 * M_PI * tStep * DT);
            IndexType NF = sourceFC.size(); 
            ComplexValueType j(0.0, 1.0); 
            if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                for (IndexType jf = 0; jf < NF; jf++) {
                    temp = scai::lama::cast<ComplexValueType>(omega);
                    temp *= -j * sourceFC[jf];
                    temp.unaryOp(temp, common::UnaryOp::EXP);
                    ftemp1 = temp;
                    ftemp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefEX());
                    fEXforward.setColumn(ftemp1, jf, common::BinaryOp::ADD);
                    ftemp1 = temp;
                    ftemp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefEY());
                    fEYforward.setColumn(ftemp1, jf, common::BinaryOp::ADD);
                }
            } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                for (IndexType jf = 0; jf < NF; jf++) {
                    temp = scai::lama::cast<ComplexValueType>(omega);
                    temp *= -j * sourceFC[jf];
                    temp.unaryOp(temp, common::UnaryOp::EXP);
                    ftemp1 = temp;
                    ftemp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefEX());
                    fEXadjoint.setColumn(ftemp1, jf, common::BinaryOp::ADD);
                    ftemp1 = temp;
                    ftemp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefEY());
                    fEYadjoint.setColumn(ftemp1, jf, common::BinaryOp::ADD);
                }
            }
        }
    } 
}

/*! \brief Sum wavefields in the frequency domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sourceFC, ValueType DT, scai::IndexType shotNumber, std::vector<scai::lama::SparseVector<ValueType>> taperEncode)
{
    double start_t_shot, end_t_shot; /* For timing */
    start_t_shot = common::Walltime::get();
    
    scai::lama::DenseMatrix<ComplexValueType> ftempM;
    scai::lama::DenseVector<ComplexValueType> ftemp;
    scai::lama::DenseVector<ComplexValueType> ftemp1;
    scai::lama::DenseVector<ComplexValueType> ftemp2;
    scai::lama::DenseVector<ValueType> xcorrSigmaEMstep;
    scai::lama::DenseVector<ValueType> xcorrEpsilonEMstep;
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
            auto dist = EXforward.getRowDistributionPtr();
            auto distNFFT = std::make_shared<scai::dmemo::NoDistribution>(nFFT);
            fEXforward = scai::lama::cast<ComplexValueType>(EXforward);
            fEXadjoint = scai::lama::cast<ComplexValueType>(EXadjoint);
            fEXforward.resize(dist, distNFFT);
            fEXadjoint.resize(dist, distNFFT);
            scai::lama::fft<ComplexValueType>(fEXforward, 1);
            scai::lama::fft<ComplexValueType>(fEXadjoint, 1);
            
            fEYforward = scai::lama::cast<ComplexValueType>(EYforward);
            fEYadjoint = scai::lama::cast<ComplexValueType>(EYadjoint);
            fEYforward.resize(dist, distNFFT);
            fEYadjoint.resize(dist, distNFFT);
            scai::lama::fft<ComplexValueType>(fEYforward, 1);
            scai::lama::fft<ComplexValueType>(fEYadjoint, 1);
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
            ftempM = scai::lama::cast<ComplexValueType>(EXforward);
            fEXforward = ftempM * L; // NR * NT * NT * NF = NR * NF
            ftempM = scai::lama::cast<ComplexValueType>(EXadjoint);
            fEXadjoint = ftempM * L;
            ftempM = scai::lama::cast<ComplexValueType>(EYforward);
            fEYforward = ftempM * L; // NR * NT * NT * NF = NR * NF
            ftempM = scai::lama::cast<ComplexValueType>(EYadjoint);
            fEYadjoint = ftempM * L;
        }            
        
        IndexType frequencySkip = 1;
        if (snapType > 0 && useSourceEncode != 0)
            frequencySkip = 2;
        
        if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fEXforward.getColumn(ftemp1, fc12Ind[jf]);
                fEXadjoint.getColumn(ftemp2, fc12Ind[jf]);                
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);
                ftemp = ftemp1;
                ftemp *= ftemp2;
                
                fEYforward.getColumn(ftemp1, fc12Ind[jf]);
                fEYadjoint.getColumn(ftemp2, fc12Ind[jf]);
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);
                ftemp1 *= ftemp2;
                ftemp1 += ftemp;
                
                xcorrSigmaEMstep = scai::lama::real(ftemp1);
                
                if (taperEncode.size() == 1) {
                    xcorrSigmaEMstep *= taperEncode[0];
                } else if (taperEncode.size() > 1) {
                    xcorrSigmaEMstep *= taperEncode[jf];
                }
                if (normalizeGradient && useSourceEncode != 0 && xcorrSigmaEMstep.maxNorm() != 0)
                    xcorrSigmaEMstep *= 1.0 / xcorrSigmaEMstep.maxNorm();
                xcorrSigmaEMstep *= weightingFreq[fcInd[jf]];
                if (snapType > 0 && workflow.workflowStage == 0 && workflow.iteration == 0) {
                    this->writeWavefield(xcorrSigmaEMstep, "xcorrSigmaEM.step", filename, frequencySkip*fcInd[jf]);
                }
                xcorrSigmaEM += xcorrSigmaEMstep;
            }
        }
        if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fEXforward.getColumn(ftemp1, fc12Ind[jf]);
                fEXadjoint.getColumn(ftemp2, fc12Ind[jf]);                
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);
                ftemp = ftemp1;
                ftemp *= ftemp2;
                
                fEYforward.getColumn(ftemp1, fc12Ind[jf]);
                fEYadjoint.getColumn(ftemp2, fc12Ind[jf]);
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);
                ftemp1 *= ftemp2;
                ftemp1 += ftemp;
                
                ftemp1 *= j * omega[jf];                
                xcorrEpsilonEMstep = scai::lama::real(ftemp1);
                
                if (taperEncode.size() == 1) {
                    xcorrEpsilonEMstep *= taperEncode[0];
                } else if (taperEncode.size() > 1) {
                    xcorrEpsilonEMstep *= taperEncode[jf];
                }
                if (normalizeGradient && useSourceEncode != 0 && xcorrEpsilonEMstep.maxNorm() != 0)
                    xcorrEpsilonEMstep *= 1.0 / xcorrEpsilonEMstep.maxNorm();
                xcorrEpsilonEMstep *= weightingFreq[fcInd[jf]];
                if (snapType > 0 && workflow.workflowStage == 0 && workflow.iteration == 0) {
                    this->writeWavefield(xcorrEpsilonEMstep, "xcorrEpsilonEM.step", filename, frequencySkip*fcInd[jf]);
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
