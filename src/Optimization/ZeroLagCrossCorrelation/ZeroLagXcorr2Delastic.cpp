#include "ZeroLagXcorr2Delastic.hpp"

using namespace scai;

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.invertForDensity)
        this->initWavefield(VSum, ctx, dist);

    if ((workflow.invertForVp) || (workflow.invertForVs) || (workflow.invertForDensity)) 
        this->initWavefield(NormalStressSum, ctx, dist);

    if ((workflow.invertForVs) || (workflow.invertForDensity)){
        this->initWavefield(ShearStress, ctx, dist);
        this->initWavefield(NormalStressDiff, ctx, dist);
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::getContextPtr()
{
    return (VSum.getContextPtr());
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
    init(ctx, dist, workflow);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::write(std::string type, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &/*workflow*/)
{
    this->writeWavefield(VSum, "VSum", type, t);
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
    if (workflow.invertForDensity)
    this->resetWavefield(VSum);

    if ((workflow.invertForVs)|| (workflow.invertForDensity)) {
    this->resetWavefield(ShearStress);
    this->resetWavefield(NormalStressDiff);
    }
    
    if ((workflow.invertForVp) || (workflow.invertForVs)|| (workflow.invertForDensity)){
    this->resetWavefield(NormalStressSum);
    }
}

/*! \brief function to update the result of the zero lag cross-correlation for per timestep 
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefields allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp1;
    lama::DenseVector<ValueType> temp2;
    if ((workflow.invertForVp) || (workflow.invertForVs) || (workflow.invertForDensity)) {
        temp1 = forwardWavefield.getRefSxx() + forwardWavefield.getRefSyy();
        temp2 = adjointWavefield.getRefSxx() + adjointWavefield.getRefSyy();
        temp1 *= temp2;
        NormalStressSum += temp1;
    }

    if ((workflow.invertForVs) || (workflow.invertForDensity)){
        temp1 = forwardWavefield.getRefSxx() - forwardWavefield.getRefSyy();
        temp2 = adjointWavefield.getRefSxx() - adjointWavefield.getRefSyy();
        temp1 *= temp2;
        NormalStressDiff += temp1;
        temp1 = forwardWavefield.getRefSxy();
        temp1 *= adjointWavefield.getRefSxy();
        ShearStress += temp1;
    }

    if (workflow.invertForDensity) {
        temp1 = forwardWavefield.getRefVX();
        temp1 *= adjointWavefield.getRefVX();
        VSum += temp1;
        temp1 = forwardWavefield.getRefVY();
        temp1 *= adjointWavefield.getRefVY();
        VSum += temp1;
    }
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::getP() const
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D elastic case.");
    return (P);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<double>;
