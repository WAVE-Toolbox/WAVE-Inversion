#include "ZeroLagXcorr2Dacoustic.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D acoustic wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::ZeroLagXcorr2Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    equationType="acoustic"; 
    numDimension=2;
    init(ctx, dist, workflow);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    type = equationType+std::to_string(numDimension)+"D";
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrRho, ctx, dist);
        if (gradientType != 0 && decomposeType != 0) {
            this->initWavefield(xcorrRhoSuRu, ctx, dist);
            this->initWavefield(xcorrRhoSdRd, ctx, dist);
            this->initWavefield(xcorrRhoSuRd, ctx, dist);
            this->initWavefield(xcorrRhoSdRu, ctx, dist);
        }
    }
    if (workflow.getInvertForVp() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrLambda, ctx, dist);
        if (gradientType != 0 && decomposeType != 0) {
            this->initWavefield(xcorrLambdaSuRu, ctx, dist);
            this->initWavefield(xcorrLambdaSdRd, ctx, dist);
            this->initWavefield(xcorrLambdaSuRd, ctx, dist);
            this->initWavefield(xcorrLambdaSdRu, ctx, dist);
        }
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getContextPtr()
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrRho, "xcorrRho", filename, t);
        if (decomposeType == 0) {
            this->writeWavefield(xcorrRhostep, "xcorrRho.step", filename, t);
        } else if (gradientType == 1 && decomposeType == 1) {
            this->writeWavefield(xcorrRhoSdRu, "xcorrRho.SdRu", filename, t);
            this->writeWavefield(xcorrRhoSuRd, "xcorrRho.SuRd", filename, t);
        } else if (gradientType == 2 && decomposeType == 1) {
            this->writeWavefield(xcorrRhoSuRu, "xcorrRho.SuRu", filename, t);
            this->writeWavefield(xcorrRhoSdRd, "xcorrRho.SdRd", filename, t);
        }
    }
    if (workflow.getInvertForVp() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrLambda, "xcorrLambda", filename, t);
        if (decomposeType == 0) {
            this->writeWavefield(xcorrLambdastep, "xcorrLambda.step", filename, t);
        } else if (gradientType == 1 && decomposeType == 1) {
            this->writeWavefield(xcorrLambdaSdRu, "xcorrLambda.SdRu", filename, t);
            this->writeWavefield(xcorrLambdaSuRd, "xcorrLambda.SuRd", filename, t);
        } else if (gradientType == 2 && decomposeType == 1) {
            this->writeWavefield(xcorrLambdaSuRu, "xcorrLambda.SuRu", filename, t);
            this->writeWavefield(xcorrLambdaSdRd, "xcorrLambda.SdRd", filename, t);
        }
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrRho);
        if (gradientType != 0 && decomposeType != 0) {
            this->resetWavefield(xcorrRhoSuRu);
            this->resetWavefield(xcorrRhoSdRd);
            this->resetWavefield(xcorrRhoSuRd);
            this->resetWavefield(xcorrRhoSdRu);
        }
    }
    if (workflow.getInvertForVp() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrLambda);
        if (gradientType != 0 && decomposeType != 0) {
            this->resetWavefield(xcorrLambdaSuRu);
            this->resetWavefield(xcorrLambdaSdRd);
            this->resetWavefield(xcorrLambdaSuRd);
            this->resetWavefield(xcorrLambdaSdRu);
        }
    }
}

/*! \brief Function to update the result of the zero lag cross-correlation per timestep 
 * 
 * The zero lag cross-correlation, \f$ X \f$, is updated with the following equations where index "forw" refers to the forward propagated wavefield and "adj" to the adjoint wavefield:
 \f{eqnarray*}
   X_{\lambda} &+=& P_{\mathrm{forw}} \cdot P_{\mathrm{adj}}  \\ 
   X_{\rho} &+=& V_{x,\mathrm{forw}} \cdot V_{x,\mathrm{adj}} + V_{y,\mathrm{forw}} \cdot V_{y,\mathrm{adj}}
 \f}
 *
 * 
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    if (workflow.getInvertForVp() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientType == 0 || decomposeType == 0) {   
            // Born kernel or FWI kernel
            xcorrLambdastep = adjointWavefield.getRefP();
            xcorrLambdastep *= forwardWavefieldDerivative.getRefP();
            xcorrLambda += xcorrLambdastep;
        } else if (gradientType == 1 && decomposeType == 1) {    
            // migration kernel using up/down-going wavefields   
            xcorrLambdastep = adjointWavefield.getRefPdown();
            xcorrLambdastep *= forwardWavefieldDerivative.getRefPdown();
            xcorrLambdaSdRu += xcorrLambdastep;           
            xcorrLambdastep = adjointWavefield.getRefPup();
            xcorrLambdastep *= forwardWavefieldDerivative.getRefPup();
            xcorrLambdaSuRd += xcorrLambdastep;  
            xcorrLambda = xcorrLambdaSdRu + xcorrLambdaSuRd;
        } else if (gradientType == 2 && decomposeType == 1) {    
            // tomographic kernel using up/down-going wavefields
            xcorrLambdastep = adjointWavefield.getRefPdown();
            xcorrLambdastep *= forwardWavefieldDerivative.getRefPup();
            xcorrLambdaSuRu += xcorrLambdastep;           
            xcorrLambdastep = adjointWavefield.getRefPup();
            xcorrLambdastep *= forwardWavefieldDerivative.getRefPdown();
            xcorrLambdaSdRd += xcorrLambdastep; 
            xcorrLambda = xcorrLambdaSuRu + xcorrLambdaSdRd;
        }
    }
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientType == 0 || decomposeType == 0) { 
            // Born kernel or FWI kernel
            xcorrRhostep = forwardWavefieldDerivative.getRefVX();
            xcorrRhostep *= adjointWavefield.getRefVX();
            xcorrRho += xcorrRhostep;
            xcorrRhostep = forwardWavefieldDerivative.getRefVY();
            xcorrRhostep *= adjointWavefield.getRefVY();
            xcorrRho += xcorrRhostep;
        } else if (gradientType == 1 && decomposeType == 1) {    
            // migration kernel using up/down-going wavefields   
            xcorrRhostep = forwardWavefieldDerivative.getRefVXdown();
            xcorrRhostep *= adjointWavefield.getRefVXdown();
            xcorrRhoSdRu += xcorrRhostep;
            xcorrRhostep = forwardWavefieldDerivative.getRefVYdown();
            xcorrRhostep *= adjointWavefield.getRefVYdown();
            xcorrRhoSdRu += xcorrRhostep;  
            xcorrRhostep = forwardWavefieldDerivative.getRefVXup();
            xcorrRhostep *= adjointWavefield.getRefVXup();
            xcorrRhoSuRd += xcorrRhostep;
            xcorrRhostep = forwardWavefieldDerivative.getRefVYup();
            xcorrRhostep *= adjointWavefield.getRefVYup();
            xcorrRhoSuRd += xcorrRhostep;  
            xcorrRho = xcorrRhoSdRu + xcorrRhoSuRd;
        } else if (gradientType == 2 && decomposeType == 1) {    
            // tomographic kernel using up/down-going wavefields
            xcorrRhostep = forwardWavefieldDerivative.getRefVXdown();
            xcorrRhostep *= adjointWavefield.getRefVXup();
            xcorrRhoSdRd += xcorrRhostep;
            xcorrRhostep = forwardWavefieldDerivative.getRefVYdown();
            xcorrRhostep *= adjointWavefield.getRefVYup();
            xcorrRhoSdRd += xcorrRhostep;
            xcorrRhostep = forwardWavefieldDerivative.getRefVXup();
            xcorrRhostep *= adjointWavefield.getRefVXdown();
            xcorrRhoSuRu += xcorrRhostep;
            xcorrRhostep = forwardWavefieldDerivative.getRefVYup();
            xcorrRhostep *= adjointWavefield.getRefVYdown();
            xcorrRhoSuRu += xcorrRhostep;
            xcorrRho = xcorrRhoSdRd + xcorrRhoSuRu;
        }
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (acoustic)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Not valid in the acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getXcorrMuA() const
{
    COMMON_THROWEXCEPTION("There is no Mu Gradient in the acoustic case.");
    return (xcorrMuA);
}

//! \brief Not valid in the acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getXcorrMuB() const
{
    COMMON_THROWEXCEPTION("There is no Mu Gradient in the acoustic case.");
    return (xcorrMuB);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getXcorrMuC() const
{
    COMMON_THROWEXCEPTION("There is no Mu Gradient in the acoustic case.");
    return (xcorrMuC);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<float>;
