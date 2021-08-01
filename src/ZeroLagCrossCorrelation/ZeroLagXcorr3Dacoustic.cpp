#include "ZeroLagXcorr3Dacoustic.hpp"

using namespace scai;



/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D acoustic wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::ZeroLagXcorr3Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    equationType="acoustic"; 
    numDimension=3;
    init(ctx, dist, workflow);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    type = equationType+std::to_string(numDimension)+"D";
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrRho, ctx, dist);
    if (workflow.getInvertForVp() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrLambda, ctx, dist);
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::getContextPtr()
{
    return (xcorrRho.getContextPtr());
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param filename file name
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrRho, "xcorrRho", filename + type, t);
    if (workflow.getInvertForVp() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrLambda, "xcorrLambda", filename + type, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrRho);
    if (workflow.getInvertForVp() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrLambda);
}

/*! \brief function to update the result of the zero lag cross-correlation for per timestep 
 *  
 * The zero lag cross-correlation, \f$ X \f$, is updated with the following equations where index "forw" refers to the forward propagated wavefield and "adj" to the adjoint wavefield:
 \f{eqnarray*}
   X_{\lambda} &+=& P_{\mathrm{forw}} \cdot P_{\mathrm{adj}}  \\ 
   X_{\rho} &+=& V_{x,\mathrm{forw}} \cdot V_{x,\mathrm{adj}} + V_{y,\mathrm{forw}} \cdot V_{y,\mathrm{adj}} + V_{z,\mathrm{forw}} \cdot V_{z,\mathrm{adj}}
 \f}
 *
 * 
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;
    if (workflow.getInvertForVp() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefieldDerivative.getRefP();
        temp *= adjointWavefield.getRefP();
        xcorrLambda += temp;
    }
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefieldDerivative.getRefVX();
        temp *= adjointWavefield.getRefVX();
        xcorrRho += temp;
        temp = forwardWavefieldDerivative.getRefVY();
        temp *= adjointWavefield.getRefVY();
        xcorrRho += temp;
        temp = forwardWavefieldDerivative.getRefVZ();
        temp *= adjointWavefield.getRefVZ();
        xcorrRho += temp;
    }
}

/*! \brief Get numDimension (3)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (acoustic)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::getXcorrMuA() const
{
    COMMON_THROWEXCEPTION("There is no Mu gradient in the acoustic case.");
    return (xcorrMuA);
}

//! \brief Not valid in the acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::getXcorrMuB() const
{
    COMMON_THROWEXCEPTION("There is no Mu gradient in the acoustic case.");
    return (xcorrMuB);
}

//! \brief Not valid in the acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::getXcorrMuC() const
{
    COMMON_THROWEXCEPTION("There is no Mu gradient in the acoustic case.");
    return (xcorrMuC);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<double>;
