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
    init(ctx, dist, workflow);
}


template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.invertForDensity)
        this->initWavefield(VSum, ctx, dist);
    if ((workflow.invertForVp) || (workflow.invertForDensity))
        this->initWavefield(P, ctx, dist);
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getContextPtr()
{
    return (VSum.getContextPtr());
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
    if (workflow.invertForDensity)
    this->writeWavefield(VSum, "VSum", type, t);
    if ((workflow.invertForVp) || (workflow.invertForDensity))
    this->writeWavefield(P, "P", type, t);
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
    if (workflow.invertForDensity)
    this->resetWavefield(VSum);
    if ((workflow.invertForVp) || (workflow.invertForDensity))
    this->resetWavefield(P);
}

/*! \brief function to update the result of the zero lag cross-correlation for per timestep 
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;
    if ((workflow.invertForVp) || (workflow.invertForDensity)) {
        temp = forwardWavefield.getRefP();
        temp *= adjointWavefield.getRefP();
        P += temp;
    }
    if (workflow.invertForDensity) {
        temp = forwardWavefield.getRefVX();
        temp *= adjointWavefield.getRefVX();
        VSum += temp;
        temp = forwardWavefield.getRefVY();
        temp *= adjointWavefield.getRefVY();
        VSum += temp;
    }
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getShearStress() const
{
    COMMON_THROWEXCEPTION("There is no ShearStress in the 2D acoustic case.");
    return (ShearStress);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getNormalStressDiff() const
{
    COMMON_THROWEXCEPTION("There is no ShearStress in the 2D acoustic case.");
    return (NormalStressDiff);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::getNormalStressSum() const
{
    COMMON_THROWEXCEPTION("There is no ShearStress in the 2D acoustic case.");
    return (NormalStressSum);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<float>;
