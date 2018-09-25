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
    if (workflow.getInvertForDensity())
        this->initWavefield(xcorrRho, ctx, dist);
    if ((workflow.getInvertForVp()) || (workflow.getInvertForDensity()))
        this->initWavefield(xcorrLambda, ctx, dist);
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getContextPtr()
{
    return (xcorrRho.getContextPtr());
}


/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::write(std::string type, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity())
        this->writeWavefield(xcorrRho, "xcorrRho", type, t);
    if ((workflow.getInvertForVp()) || (workflow.getInvertForDensity()))
        this->writeWavefield(xcorrLambda, "xcorrLambda", type, t);
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::writeSnapshot(IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    write(type, t, workflow);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity())
        this->resetWavefield(xcorrRho);
    if ((workflow.getInvertForVp()) || (workflow.getInvertForDensity()))
        this->resetWavefield(xcorrLambda);
}

/*! \brief Function to update the result of the zero lag cross-correlation per timestep 
 * 
 * The zero lag cross-correlation, \f$ X \f$, is updated with the following equations where index "forw" refers to the forward propagated wavefield and "adj" to the adjoint wavefield:
 \f{eqnarray*}
   X_{\lambda} &+=& P_{\mathrm{forw}} \cdot P_{\mathrm{adj}}  \\ 
   X_{\rho} &+=& V_{x,\mathrm{forw}} \cdot V_{x,\mathrm{adj}} + V_{y,\mathrm{forw}} \cdot V_{y,\mathrm{adj}}
 \f}
 *
 *  Note that the forwardWavefield is actually the derivative of the forward wavefield (see variable wavefieldrecord in IFOS.cpp).
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;
    if ((workflow.getInvertForVp()) || (workflow.getInvertForDensity())) {
        temp = forwardWavefield.getRefP();
        temp *= adjointWavefield.getRefP();
        xcorrLambda += temp;
    }
    if (workflow.getInvertForDensity()) {
        temp = forwardWavefield.getRefVX();
        temp *= adjointWavefield.getRefVX();
        xcorrRho += temp;
        temp = forwardWavefield.getRefVY();
        temp *= adjointWavefield.getRefVY();
        xcorrRho += temp;
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
