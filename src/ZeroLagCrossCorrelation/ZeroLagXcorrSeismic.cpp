#include "ZeroLagXcorrSeismic.hpp"

using namespace scai;

//! \brief Getter routine for VSum wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrRho() const
{
    return (xcorrRho);
}

//! \brief Getter routine for XcorrLambda
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrLambda() const
{
    return (xcorrLambda);
}

//! \brief Getter routine for xcorrMuA
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrMuA() const
{
    return (xcorrMuA);
}

//! \brief Getter routine for xcorrMuB
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrMuB() const
{
    return (xcorrMuB);
}

//! \brief Getter routine for xcorrMuC
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrMuC() const
{
    return (xcorrMuC);
}


/* EM */
//! \brief Getter routine for xcorrSigma
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrSigma() const
{
    COMMON_THROWEXCEPTION("There is no xcorrSigma in the ZeroLagXcorrSeismic.")
    return (xcorrSigma);
}

//! \brief Getter routine for xcorrEpsilon
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrEpsilon() const
{
    COMMON_THROWEXCEPTION("There is no xcorrEpsilon in the ZeroLagXcorrSeismic.")
    return (xcorrEpsilon);
}

//! \brief Getter routine for xcorrREpsilonSigma
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrREpsilonSigma() const
{
    COMMON_THROWEXCEPTION("There is no xcorrREpsilonSigma in the ZeroLagXcorrSeismic.")
    return (xcorrREpsilonSigma);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<double>;
