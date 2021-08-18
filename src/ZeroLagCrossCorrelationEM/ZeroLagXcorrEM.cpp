#include "ZeroLagXcorrEM.hpp"

using namespace scai;

//! \brief Getter routine for xcorrSigmaEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrSigmaEM() const
{
    return (xcorrSigmaEM);
}

//! \brief Getter routine for xcorrEpsilonEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrEpsilonEM() const
{
    return (xcorrEpsilonEM);
}

//! \brief Getter routine for xcorrRSigmaEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrRSigmaEM() const
{
    return (xcorrRSigmaEM);
}

//! \brief Getter routine for xcorrREpsilonEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrREpsilonEM() const
{
    return (xcorrREpsilonEM);
}


/* Seismic */
//! \brief Getter routine for VSum wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrRho() const
{
    COMMON_THROWEXCEPTION("There is no xcorrRho in the ZeroLagXcorrEM.")
    return (xcorrRho);
}

//! \brief Getter routine for XcorrLambda
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrLambda() const
{
    COMMON_THROWEXCEPTION("There is no xcorrLambda in the ZeroLagXcorrEM.")
    return (xcorrLambda);
}

//! \brief Getter routine for xcorrMuA
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrMuA() const
{
    COMMON_THROWEXCEPTION("There is no xcorrMuA in the ZeroLagXcorrEM.")
    return (xcorrMuA);
}

//! \brief Getter routine for xcorrMuB
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrMuB() const
{
    COMMON_THROWEXCEPTION("There is no xcorrMuB in the ZeroLagXcorrEM.")
    return (xcorrMuB);
}

//! \brief Getter routine for xcorrMuC
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrMuC() const
{
    COMMON_THROWEXCEPTION("There is no xcorrMuC in the ZeroLagXcorrEM.")
    return (xcorrMuC);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<double>;
