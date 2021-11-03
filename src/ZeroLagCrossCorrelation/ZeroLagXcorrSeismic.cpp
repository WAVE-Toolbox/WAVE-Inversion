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
//! \brief Getter routine for xcorrSigmaEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrSigmaEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrSigmaEM in the ZeroLagXcorrSeismic.")
    return (xcorrSigmaEM);
}

//! \brief Getter routine for xcorrEpsilonEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrEpsilonEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrEpsilonEM in the ZeroLagXcorrSeismic.")
    return (xcorrEpsilonEM);
}

//! \brief Getter routine for xcorrRSigmaEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrRSigmaEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrRSigmaEM in the ZeroLagXcorrSeismic.")
    return (xcorrRSigmaEM);
}

//! \brief Getter routine for xcorrREpsilonEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<ValueType>::getXcorrREpsilonEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrREpsilonEM in the ZeroLagXcorrSeismic.")
    return (xcorrREpsilonEM);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorrSeismic<double>;