#include "ZeroLagXcorr2Delastic.hpp"

using namespace scai;

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
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::ZeroLagXcorr2Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(VSum, ctx, dist);
    COMMON_THROWEXCEPTION("elastic convolution is not implemented yet.")
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::write(std::string type, IndexType t)
{
    this->writeWavefield(VSum, "VSum", type, t);
    COMMON_THROWEXCEPTION("elastic convolution is not implemented yet.")
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::writeSnapshot(IndexType t)
{
    write(type, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::reset()
{
    this->resetWavefield(VSum);
    COMMON_THROWEXCEPTION("elastic convolution is not implemented yet.")
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::getP() const
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D elastic case.")
    return (P);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<double>;
