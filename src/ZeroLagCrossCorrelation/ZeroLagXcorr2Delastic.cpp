#include "ZeroLagXcorr2Delastic.hpp"

using namespace scai;

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity())
        this->initWavefield(xcorrRho, ctx, dist);

    if ((workflow.getInvertForVp()) || (workflow.getInvertForVs()) || (workflow.getInvertForDensity()))
        this->initWavefield(xcorrLambda, ctx, dist);

    if ((workflow.getInvertForVs()) || (workflow.getInvertForDensity())) {
        this->initWavefield(xcorrMuA, ctx, dist);
        this->initWavefield(xcorrMuB, ctx, dist);
        this->initWavefield(xcorrMuC, ctx, dist);
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::getContextPtr()
{
    return (xcorrRho.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D elastic wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::ZeroLagXcorr2Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    equationType="elastic"; 
    numDimension=2;
    init(ctx, dist, workflow);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::write(std::string type, IndexType t, KITGPI::Workflow::Workflow<ValueType> const & /*workflow*/)
{
    this->writeWavefield(xcorrRho, "xcorrRho", type, t);
    COMMON_THROWEXCEPTION("write correlated elastic wavefields is not implemented yet.")
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::writeSnapshot(IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    write(type, t, workflow);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForDensity())
        this->resetWavefield(xcorrRho);

    if ((workflow.getInvertForVs()) || (workflow.getInvertForDensity())) {
        this->resetWavefield(xcorrMuA);
        this->resetWavefield(xcorrMuB);
        this->resetWavefield(xcorrMuC);
    }

    if ((workflow.getInvertForVp()) || (workflow.getInvertForVs()) || (workflow.getInvertForDensity())) {
        this->resetWavefield(xcorrLambda);
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
 * Note that the forwardWavefield is actually the derivative of the forward wavefield (see variable wavefieldrecord in main.cpp).
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefields allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp1;
    lama::DenseVector<ValueType> temp2;
    if ((workflow.getInvertForVp()) || (workflow.getInvertForVs()) || (workflow.getInvertForDensity())) {
        temp1 = forwardWavefield.getRefSxx() + forwardWavefield.getRefSyy();
        temp2 = adjointWavefield.getRefSxx() + adjointWavefield.getRefSyy();
        temp1 *= temp2;
        xcorrLambda += temp1;
    }

    if ((workflow.getInvertForVs()) || (workflow.getInvertForDensity())) {
        temp1 = forwardWavefield.getRefSxx();
        temp1 *= adjointWavefield.getRefSxx();
        xcorrMuA += temp1;
        temp1 = forwardWavefield.getRefSyy();
        temp1 *= adjointWavefield.getRefSyy();
        xcorrMuA += temp1;

        temp1 = forwardWavefield.getRefSyy();
        temp1 *= adjointWavefield.getRefSxx();
        xcorrMuB += temp1;
        temp1 = forwardWavefield.getRefSxx();
        temp1 *= adjointWavefield.getRefSyy();
        xcorrMuB += temp1;

        temp1 = forwardWavefield.getRefSxy();
        temp1 *= adjointWavefield.getRefSxy();
        xcorrMuC += temp1;
    }

    if (workflow.getInvertForDensity()) {
        temp1 = forwardWavefield.getRefVX();
        temp1 *= adjointWavefield.getRefVX();
        xcorrRho += temp1;
        temp1 = forwardWavefield.getRefVY();
        temp1 *= adjointWavefield.getRefVY();
        xcorrRho += temp1;
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
