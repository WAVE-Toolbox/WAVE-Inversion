#include "ZeroLagXcorr2Delastic.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D elastic wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::ZeroLagXcorr2Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot)
{
    equationType="elastic"; 
    numDimension=2;
    init(ctx, dist, workflow, config, numShotPerSuperShot);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot)
{
    type = equationType+std::to_string(numDimension)+"D";
    gradientKernel = config.getAndCatch("gradientKernel", 0);
    decomposition = config.getAndCatch("decomposition", 0);
    useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    gradientDomain = config.getAndCatch("gradientDomain", 0);
    
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrRho, ctx, dist);
    }

    if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrLambda, ctx, dist);
    }

    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrMuA, ctx, dist);
        this->initWavefield(xcorrMuB, ctx, dist);
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
                VXforward.setContextPtr(ctx);
                VXadjoint.setContextPtr(ctx);
                VXforward.allocate(dist, no_dist_NT);
                VXadjoint.allocate(dist, no_dist_NT);
                VYforward.setContextPtr(ctx);
                VYadjoint.setContextPtr(ctx);
                VYforward.allocate(dist, no_dist_NT);
                VYadjoint.allocate(dist, no_dist_NT);
            }
            if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                Sxxforward.setContextPtr(ctx);
                Sxxadjoint.setContextPtr(ctx);
                Sxxforward.allocate(dist, no_dist_NT);
                Sxxadjoint.allocate(dist, no_dist_NT);
                Syyforward.setContextPtr(ctx);
                Syyadjoint.setContextPtr(ctx);
                Syyforward.allocate(dist, no_dist_NT);
                Syyadjoint.allocate(dist, no_dist_NT);
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                Sxyforward.setContextPtr(ctx);
                Sxyadjoint.setContextPtr(ctx);
                Sxyforward.allocate(dist, no_dist_NT);
                Sxyadjoint.allocate(dist, no_dist_NT);
            }
        } else if (gradientDomain == 3) {
            dmemo::DistributionPtr no_dist_NS(new scai::dmemo::NoDistribution(numShotPerSuperShot));
            if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                fVXforward.setContextPtr(ctx);
                fVXadjoint.setContextPtr(ctx);
                fVXforward.allocate(dist, no_dist_NS);
                fVXadjoint.allocate(dist, no_dist_NS);
                fVYforward.setContextPtr(ctx);
                fVYadjoint.setContextPtr(ctx);
                fVYforward.allocate(dist, no_dist_NS);
                fVYadjoint.allocate(dist, no_dist_NS);
            }
            if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                fSxxforward.setContextPtr(ctx);
                fSxxadjoint.setContextPtr(ctx);
                fSxxforward.allocate(dist, no_dist_NS);
                fSxxadjoint.allocate(dist, no_dist_NS);
                fSyyforward.setContextPtr(ctx);
                fSyyadjoint.setContextPtr(ctx);
                fSyyforward.allocate(dist, no_dist_NS);
                fSyyadjoint.allocate(dist, no_dist_NS);
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                fSxyforward.setContextPtr(ctx);
                fSxyadjoint.setContextPtr(ctx);
                fSxyforward.allocate(dist, no_dist_NS);
                fSxyadjoint.allocate(dist, no_dist_NS);
            }
        }
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::getContextPtr()
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrRho, "xcorrRho", filename, t);
    }

    if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrLambda, "xcorrLambda", filename, t);
    }

    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrMuA, "xcorrMuA", filename, t);
        this->writeWavefield(xcorrMuB, "xcorrMuB", filename, t);
        this->writeWavefield(xcorrMuC, "xcorrMuC", filename, t);
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrRho);
        if (gradientDomain == 1 || gradientDomain == 2) {
            VXforward.scale(0.0);
            VXadjoint.scale(0.0);
            VYforward.scale(0.0);
            VYadjoint.scale(0.0);
        } else if (gradientDomain == 3) {
            fVXforward.scale(0.0);
            fVXadjoint.scale(0.0);
            fVYforward.scale(0.0);
            fVYadjoint.scale(0.0);
        }
    }

    if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrLambda);
        if (gradientDomain == 1 || gradientDomain == 2) {
            Sxxforward.scale(0.0);
            Sxxadjoint.scale(0.0);
            Syyforward.scale(0.0);
            Syyadjoint.scale(0.0);
        } else if (gradientDomain == 3) {
            fSxxforward.scale(0.0);
            fSxxadjoint.scale(0.0);
            fSyyforward.scale(0.0);
            fSyyadjoint.scale(0.0);
        }
    }

    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrMuA);
        this->resetWavefield(xcorrMuB);
        this->resetWavefield(xcorrMuC);
        if (gradientDomain == 1 || gradientDomain == 2) {
            Sxyforward.scale(0.0);
            Sxyadjoint.scale(0.0);
        } else if (gradientDomain == 3) {
            fSxyforward.scale(0.0);
            fSxyadjoint.scale(0.0);
        }
    }
}

/*! \brief function to update the result of the zero lag cross-correlation for per timestep 
 * 
 * The zero lag cross-correlation, \f$ X \f$, is updated with the following equations where index "forw" refers to the forward propagated wavefield and "adj" to the adjoint wavefield. \f$ S_{ii} \f$ are the stresses and \f$ V_i \f$ the particle velocities.
 \f{eqnarray*}
   X_{\lambda} &+=& (S_{xx,\mathrm{forw}} + S_{yy,\mathrm{forw}}) \cdot (S_{xx,\mathrm{adj}} + S_{yy,\mathrm{adj}})  \\ 
   X_{\mu,A} &+=& S_{xx,\mathrm{forw}} \cdot S_{xx,\mathrm{adj}} + S_{yy,\mathrm{forw}} \cdot S_{yy,\mathrm{adj}} \\
   X_{\mu,B} &+=& S_{yy,\mathrm{forw}} \cdot S_{xx,\mathrm{adj}} + S_{xx,\mathrm{forw}} \cdot S_{yy,\mathrm{adj}} \\
   X_{\mu,C} &+=& S_{xy,\mathrm{forw}} \cdot S_{xy,\mathrm{adj}}  \\
   
   X_{\rho} &+=& V_{x,\mathrm{forw}} \cdot V_{x,\mathrm{adj}} + V_{y,\mathrm{forw}} \cdot V_{y,\mathrm{adj}}
 \f}
 * 
 * 
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> temp;
    scai::lama::DenseVector<ValueType> temp1;
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefieldDerivative.getRefVX();
        temp *= adjointWavefield.getRefVX();
        xcorrRho += temp;
        temp = forwardWavefieldDerivative.getRefVY();
        temp *= adjointWavefield.getRefVY();
        xcorrRho += temp;
    }
    
    if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefieldDerivative.getRefSxx() + forwardWavefieldDerivative.getRefSyy();
        temp1 = adjointWavefield.getRefSxx() + adjointWavefield.getRefSyy();
        temp *= temp1;
        xcorrLambda += temp;
    }

    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefieldDerivative.getRefSxx();
        temp *= adjointWavefield.getRefSxx();
        xcorrMuA += temp;
        temp = forwardWavefieldDerivative.getRefSyy();
        temp *= adjointWavefield.getRefSyy();
        xcorrMuA += temp;

        temp = forwardWavefieldDerivative.getRefSyy();
        temp *= adjointWavefield.getRefSxx();
        xcorrMuB += temp;
        temp = forwardWavefieldDerivative.getRefSxx();
        temp *= adjointWavefield.getRefSyy();
        xcorrMuB += temp;

        temp = forwardWavefieldDerivative.getRefSxy();
        temp *= adjointWavefield.getRefSxy();
        xcorrMuC += temp;
    }
}

/*! \brief Gather wavefields in the time domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::gatherWavefields(Wavefields::Wavefields<ValueType> &wavefields, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint)
{
    IndexType T0 = 0;
    if (isAdjoint || useSourceEncode != 0)
        T0 = NT; // half recording time.
    if (gradientKernel == 0 || decomposition == 0) {   
        if (gradientDomain == 1 || gradientDomain == 2) {
            if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                    VXforward.setColumn(wavefields.getRefVX(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
                    VYforward.setColumn(wavefields.getRefVY(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
                } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                    VXadjoint.setColumn(wavefields.getRefVX(), tStep / workflow.skipDT, common::BinaryOp::COPY);
                    VYadjoint.setColumn(wavefields.getRefVY(), tStep / workflow.skipDT, common::BinaryOp::COPY);
                }
            }
            if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                    Sxxforward.setColumn(wavefields.getRefSxx(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
                    Syyforward.setColumn(wavefields.getRefSyy(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
                } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                    Sxxadjoint.setColumn(wavefields.getRefSxx(), tStep / workflow.skipDT, common::BinaryOp::COPY);
                    Syyadjoint.setColumn(wavefields.getRefSyy(), tStep / workflow.skipDT, common::BinaryOp::COPY);
                }
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                    Sxyforward.setColumn(wavefields.getRefSxy(), tStep / workflow.skipDT - T0, common::BinaryOp::COPY);
                } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                    Sxyadjoint.setColumn(wavefields.getRefSxy(), tStep / workflow.skipDT, common::BinaryOp::COPY);
                }
            }
        } else if (gradientDomain == 3) {
            scai::lama::DenseVector<ComplexValueType> temp;
            scai::lama::DenseVector<ComplexValueType> temp1;
            scai::lama::DenseVector<ValueType> omega(wavefields.getRefVX().getDistributionPtr(), 2.0 * M_PI * tStep * DT);
            IndexType NF = sourceFC.size(); 
            ComplexValueType j(0.0, 1.0); 
            if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                    for (IndexType jf = 0; jf < NF; jf++) {
                        temp = scai::lama::cast<ComplexValueType>(omega);
                        temp *= -j * sourceFC[jf];
                        temp.unaryOp(temp, common::UnaryOp::EXP);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefVX());
                        fVXforward.setColumn(temp1, jf, common::BinaryOp::ADD);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefVY());
                        fVYforward.setColumn(temp1, jf, common::BinaryOp::ADD);
                    }
                } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                    for (IndexType jf = 0; jf < NF; jf++) {
                        temp = scai::lama::cast<ComplexValueType>(omega);
                        temp *= -j * sourceFC[jf];
                        temp.unaryOp(temp, common::UnaryOp::EXP);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefVX());
                        fVXadjoint.setColumn(temp1, jf, common::BinaryOp::ADD);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefVY());
                        fVYadjoint.setColumn(temp1, jf, common::BinaryOp::ADD);
                    }
                }
            }
            if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                    for (IndexType jf = 0; jf < NF; jf++) {
                        temp = scai::lama::cast<ComplexValueType>(omega);
                        temp *= -j * sourceFC[jf];
                        temp.unaryOp(temp, common::UnaryOp::EXP);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefSxx());
                        fSxxforward.setColumn(temp1, jf, common::BinaryOp::ADD);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefSyy());
                        fSyyforward.setColumn(temp1, jf, common::BinaryOp::ADD);
                    }
                } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                    for (IndexType jf = 0; jf < NF; jf++) {
                        temp = scai::lama::cast<ComplexValueType>(omega);
                        temp *= -j * sourceFC[jf];
                        temp.unaryOp(temp, common::UnaryOp::EXP);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefSxx());
                        fSxxadjoint.setColumn(temp1, jf, common::BinaryOp::ADD);
                        temp1 = temp;
                        temp1 *= scai::lama::cast<ComplexValueType>(wavefields.getRefSyy());
                        fSyyadjoint.setColumn(temp1, jf, common::BinaryOp::ADD);
                    }
                }
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                if (!isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT >= T0))) {
                    for (IndexType jf = 0; jf < NF; jf++) {
                        temp = scai::lama::cast<ComplexValueType>(omega);
                        temp *= -j * sourceFC[jf];
                        temp.unaryOp(temp, common::UnaryOp::EXP);
                        temp *= scai::lama::cast<ComplexValueType>(wavefields.getRefSxy());
                        fSxyforward.setColumn(temp, jf, common::BinaryOp::ADD);
                    }
                } else if (isAdjoint && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep / workflow.skipDT < T0))) {
                    for (IndexType jf = 0; jf < NF; jf++) {
                        temp = scai::lama::cast<ComplexValueType>(omega);
                        temp *= -j * sourceFC[jf];
                        temp.unaryOp(temp, common::UnaryOp::EXP);
                        temp *= scai::lama::cast<ComplexValueType>(wavefields.getRefSxy());
                        fSxyforward.setColumn(temp, jf, common::BinaryOp::ADD);
                    }
                }
            }
        }
    } 
}

/*! \brief Sum wavefields in the frequency domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sourceFC, ValueType DT, scai::IndexType shotNumber, std::vector<scai::lama::SparseVector<ValueType>> taperEncode)
{
    double start_t_shot, end_t_shot; /* For timing */
    start_t_shot = common::Walltime::get();
    
    scai::lama::DenseMatrix<ComplexValueType> ftempM;
    scai::lama::DenseVector<ComplexValueType> ftemp;
    scai::lama::DenseVector<ComplexValueType> ftemp1;
    scai::lama::DenseVector<ComplexValueType> ftemp2;
    scai::lama::DenseVector<ValueType> xcorrRhostep;
    scai::lama::DenseVector<ValueType> xcorrLambdastep;
    scai::lama::DenseVector<ValueType> xcorrMuAstep;
    scai::lama::DenseVector<ValueType> xcorrMuBstep;
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
            auto dist = VXforward.getRowDistributionPtr();
            auto distNFFT = std::make_shared<scai::dmemo::NoDistribution>(nFFT);
            if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                fVXforward = scai::lama::cast<ComplexValueType>(VXforward);
                fVXadjoint = scai::lama::cast<ComplexValueType>(VXadjoint);
                fVXforward.resize(dist, distNFFT);
                fVXadjoint.resize(dist, distNFFT);
                scai::lama::fft<ComplexValueType>(fVXforward, 1);
                scai::lama::fft<ComplexValueType>(fVXadjoint, 1);

                fVYforward = scai::lama::cast<ComplexValueType>(VYforward);
                fVYadjoint = scai::lama::cast<ComplexValueType>(VYadjoint);
                fVYforward.resize(dist, distNFFT);
                fVYadjoint.resize(dist, distNFFT);
                scai::lama::fft<ComplexValueType>(fVYforward, 1);
                scai::lama::fft<ComplexValueType>(fVYadjoint, 1);
            }
            if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                fSxxforward = scai::lama::cast<ComplexValueType>(Sxxforward);
                fSxxadjoint = scai::lama::cast<ComplexValueType>(Sxxadjoint);
                fSxxforward.resize(dist, distNFFT);
                fSxxadjoint.resize(dist, distNFFT);
                scai::lama::fft<ComplexValueType>(fSxxforward, 1);
                scai::lama::fft<ComplexValueType>(fSxxadjoint, 1);

                fSyyforward = scai::lama::cast<ComplexValueType>(Syyforward);
                fSyyadjoint = scai::lama::cast<ComplexValueType>(Syyadjoint);
                fSyyforward.resize(dist, distNFFT);
                fSyyadjoint.resize(dist, distNFFT);
                scai::lama::fft<ComplexValueType>(fSyyforward, 1);
                scai::lama::fft<ComplexValueType>(fSyyadjoint, 1);
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                fSxyforward = scai::lama::cast<ComplexValueType>(Sxyforward);
                fSxyadjoint = scai::lama::cast<ComplexValueType>(Sxyadjoint);
                fSxyforward.resize(dist, distNFFT);
                fSxyadjoint.resize(dist, distNFFT);
                scai::lama::fft<ComplexValueType>(fSxyforward, 1);
                scai::lama::fft<ComplexValueType>(fSxyadjoint, 1);
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
                ftempM = scai::lama::cast<ComplexValueType>(VXforward);
                fVXforward = ftempM * L; // NR * NT * NT * NF = NR * NF
                ftempM = scai::lama::cast<ComplexValueType>(VXadjoint);
                fVXadjoint = ftempM * L;
                ftempM = scai::lama::cast<ComplexValueType>(VYforward);
                fVYforward = ftempM * L; // NR * NT * NT * NF = NR * NF
                ftempM = scai::lama::cast<ComplexValueType>(VYadjoint);
                fVYadjoint = ftempM * L;
            }
            if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                ftempM = scai::lama::cast<ComplexValueType>(Sxxforward);
                fSxxforward = ftempM * L; // NR * NT * NT * NF = NR * NF
                ftempM = scai::lama::cast<ComplexValueType>(Sxxadjoint);
                fSxxadjoint = ftempM * L;
                ftempM = scai::lama::cast<ComplexValueType>(Syyforward);
                fSyyforward = ftempM * L; // NR * NT * NT * NF = NR * NF
                ftempM = scai::lama::cast<ComplexValueType>(Syyadjoint);
                fSyyadjoint = ftempM * L;
            }
            if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
                ftempM = scai::lama::cast<ComplexValueType>(Sxyforward);
                fSxyforward = ftempM * L; // NR * NT * NT * NF = NR * NF
                ftempM = scai::lama::cast<ComplexValueType>(Sxyadjoint);
                fSxyadjoint = ftempM * L;
            }
        }            
        
        IndexType frequencySkip = 1;
        if (snapType > 0 && useSourceEncode != 0)
            frequencySkip = 2;
        
        if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {                
                fVXforward.getColumn(ftemp1, fc12Ind[jf]);
                fVXadjoint.getColumn(ftemp2, fc12Ind[jf]);              
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);                  
                ftemp = ftemp1;               
                ftemp *= ftemp2;
                            
                fVYforward.getColumn(ftemp1, fc12Ind[jf]);
                fVYadjoint.getColumn(ftemp2, fc12Ind[jf]);              
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);                
                ftemp1 *= ftemp2;
                ftemp1 += ftemp;
                
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
        if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fSxxforward.getColumn(ftemp1, fc12Ind[jf]);  
                fSyyforward.getColumn(ftemp2, fc12Ind[jf]);                     
                ftemp = ftemp1;
                ftemp += ftemp2;
                
                fSxxadjoint.getColumn(ftemp1, fc12Ind[jf]);  
                fSyyadjoint.getColumn(ftemp2, fc12Ind[jf]); 
                ftemp1 += ftemp2;
                ftemp1.unaryOp(ftemp1, common::UnaryOp::CONJ);  
                ftemp1 *= ftemp;
                
                xcorrLambdastep = scai::lama::real(ftemp1);
                
                if (normalizeGradient && useSourceEncode != 0 && xcorrLambdastep.maxNorm() != 0)
                    xcorrLambdastep *= 1.0 / xcorrLambdastep.maxNorm();
                xcorrLambdastep *= weightingFreq[fcInd[jf]];
                if (snapType > 0 && workflow.workflowStage == 0 && workflow.iteration == 0) {
                    this->writeWavefield(xcorrLambdastep, "xcorrLambda.step", filename, frequencySkip*fcInd[jf]);
                }
                xcorrLambda += xcorrLambdastep;
            }
        }
        if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
            for (IndexType jf = 0; jf < nfc12; jf++) {
                fSxxforward.getColumn(ftemp1, fc12Ind[jf]);
                fSxxadjoint.getColumn(ftemp2, fc12Ind[jf]);               
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);                
                ftemp = ftemp1;
                ftemp *= ftemp2;
                
                fSyyforward.getColumn(ftemp1, fc12Ind[jf]);
                fSyyadjoint.getColumn(ftemp2, fc12Ind[jf]); 
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);
                ftemp1 *= ftemp2;
                ftemp1 += ftemp;
                
                xcorrMuAstep = scai::lama::real(ftemp1);
                
                fSxxforward.getColumn(ftemp1, fc12Ind[jf]);
                fSyyadjoint.getColumn(ftemp2, fc12Ind[jf]);               
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);                
                ftemp = ftemp1;
                ftemp *= ftemp2;
                
                fSyyforward.getColumn(ftemp1, fc12Ind[jf]);
                fSxxadjoint.getColumn(ftemp2, fc12Ind[jf]); 
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);
                ftemp1 *= ftemp2;
                ftemp1 += ftemp;
                
                xcorrMuBstep = scai::lama::real(ftemp1);
                
                fSxyforward.getColumn(ftemp1, fc12Ind[jf]);
                fSxyadjoint.getColumn(ftemp2, fc12Ind[jf]); 
                ftemp2.unaryOp(ftemp2, common::UnaryOp::CONJ);
                ftemp1 *= ftemp2;
                
                xcorrMuCstep = scai::lama::real(ftemp1);
                
                if (normalizeGradient && useSourceEncode != 0 && xcorrMuCstep.maxNorm() != 0) {
                    xcorrMuAstep *= 1.0 / xcorrMuAstep.maxNorm();
                    xcorrMuBstep *= 1.0 / xcorrMuBstep.maxNorm();
                    xcorrMuCstep *= 1.0 / xcorrMuCstep.maxNorm();
                }
                xcorrMuAstep *= weightingFreq[fcInd[jf]];
                xcorrMuBstep *= weightingFreq[fcInd[jf]];
                xcorrMuCstep *= weightingFreq[fcInd[jf]];
                if (snapType > 0 && workflow.workflowStage == 0 && workflow.iteration == 0) {
                    this->writeWavefield(xcorrMuAstep, "xcorrMuA.step", filename, frequencySkip*fcInd[jf]);
                    this->writeWavefield(xcorrMuBstep, "xcorrMuB.step", filename, frequencySkip*fcInd[jf]);
                    this->writeWavefield(xcorrMuCstep, "xcorrMuC.step", filename, frequencySkip*fcInd[jf]);
                }
                xcorrMuA += xcorrMuAstep;
                xcorrMuB += xcorrMuBstep;
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {       
        xcorrRho = lhs * xcorrRho;
    }
    if (workflow.getInvertForVp() || workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        xcorrLambda = lhs * xcorrLambda;
    }
    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        xcorrMuA = lhs * xcorrMuA;
        xcorrMuB = lhs * xcorrMuB;
        xcorrMuC = lhs * xcorrMuC;
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (elastic)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::getEquationType() const
{
    return (equationType);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<double>;
