#include "ZeroLagXcorr2Dviscosh.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D viscosh wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::ZeroLagXcorr2Dviscosh(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    equationType="viscosh"; 
    numDimension=2;
    init(ctx, dist, workflow);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    type = equationType+std::to_string(numDimension)+"D";
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrRho, ctx, dist);
    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrMuC, ctx, dist);
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::getContextPtr()
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrRho, "xcorrRho", filename, t);
    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrMuC, "xcorrMuC", filename, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrRho);
    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrMuC);
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;
    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefieldDerivative.getRefSxz();
        temp *= adjointWavefield.getRefSxz();
        xcorrMuC += temp;
        temp = forwardWavefieldDerivative.getRefSyz();
        temp *= adjointWavefield.getRefSyz();
        xcorrMuC += temp;
    }
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefieldDerivative.getRefVZ();
        temp *= adjointWavefield.getRefVZ();
        xcorrRho += temp;
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (viscosh)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Not valid in the viscosh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::getXcorrMuA() const
{
    COMMON_THROWEXCEPTION("There is no MuA Gradient in the viscosh case.");
    return (xcorrMuA);
}

//! \brief Not valid in the viscosh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::getXcorrMuB() const
{
    COMMON_THROWEXCEPTION("There is no MuB Gradient in the viscosh case.");
    return (xcorrMuB);
}

//! \brief Not valid in the viscosh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<ValueType>::getXcorrLambda() const
{
    COMMON_THROWEXCEPTION("There is no Lambda Gradient in the viscosh case.");
    return (xcorrLambda);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscosh<float>;