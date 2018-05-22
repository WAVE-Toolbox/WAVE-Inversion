#include "ZeroLagXcorr3Delastic.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<ValueType>::getContextPtr()
{
    return (VSum.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D elastic wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<ValueType>::ZeroLagXcorr3Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    init(ctx, dist, workflow);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &/*workflow*/)
{
    this->initWavefield(VSum, ctx, dist);
    COMMON_THROWEXCEPTION("3Delastic convolution is not implemented yet.")
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<ValueType>::write(std::string type, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &/*workflow*/)
{
    this->writeWavefield(VSum, "VSum", type, t);
    COMMON_THROWEXCEPTION("3Delastic convolution is not implemented yet.")
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<ValueType>::writeSnapshot(IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    write(type, t, workflow);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &/*workflow*/)
{
    this->resetWavefield(VSum);
    COMMON_THROWEXCEPTION("3Delastic convolution is not implemented yet.")
}

//! \brief Not valid in the 3D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<ValueType>::getP() const
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 3D elastic case.")
    return (P);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<double>;
