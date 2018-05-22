#include "ZeroLagXcorr3Dacoustic.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::getContextPtr()
{
    return (VSum.getContextPtr());
}

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
    init(ctx, dist, workflow);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &/*workflow*/)
{
    this->initWavefield(VSum, ctx, dist);
    this->initWavefield(P, ctx, dist);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::write(std::string type, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &/*workflow*/)
{
    this->writeWavefield(VSum, "VSum", type, t);
    this->writeWavefield(P, "P", type, t);
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::writeSnapshot(IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    write(type, t, workflow);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &/*workflow*/)
{
    this->resetWavefield(VSum);
    this->resetWavefield(P);
}

/*! \brief function to update the result of the zero lag cross-correlation for per timestep 
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::update(Wavefields::Wavefields<ValueType> &/*forwardWavefield*/, Wavefields::Wavefields<ValueType> &/*adjointWavefield*/, KITGPI::Workflow::Workflow<ValueType> const &/*workflow*/)
{
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<double>;
