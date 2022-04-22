#include "ZeroLagXcorrEM.hpp"

using namespace scai;

//! \brief Getter routine for xcorrSigma
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrSigma() const
{
    return (xcorrSigma);
}

//! \brief Getter routine for xcorrEpsilon
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrEpsilon() const
{
    return (xcorrEpsilon);
}

//! \brief Getter routine for xcorrREpsilonSigma
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::getXcorrREpsilonSigma() const
{
    return (xcorrREpsilonSigma);
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
