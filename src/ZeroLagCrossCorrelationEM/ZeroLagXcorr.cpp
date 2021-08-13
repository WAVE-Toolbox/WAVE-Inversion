#include "ZeroLagXcorr.hpp"

using namespace scai;

/*! \brief Reset a single wavefield to zero.
 \param vector Vector to be reset to 0
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::resetWavefield(scai::lama::DenseVector<ValueType> &vector)
{
    vector = 0;
}

/*! \brief set decomposeType
 \param setDecomposeType setDecomposeType
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::setDecomposeType(scai::IndexType setDecomposeType)
{
    decomposeType = setDecomposeType;
}

/*! \brief set gradientType
 \param setGradientType setGradientType
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::setGradientType(scai::IndexType setGradientType)
{
    gradientType = setGradientType;
}

/*! \brief Intitialisation of a single wavefield vector.
 *
 * This method will set the context, allocate the the wavefield and set the field to zero.
 *
 \param vector Vector to be set
 \param ctx Context pointer
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    // NOTE_TB: why is this an own method,

    vector = lama::zero<lama::DenseVector<ValueType>>(dist, ctx);

    // vector.setContextPtr(ctx);
    // vector.allocate(dist);
    // resetWavefield(vector);
}

/*! \brief Methode to Write Wavefield for timestep t
 *
 \param vector Vector written to file
 \param type Wavefield-type ( elastic, viscoelastic)
 \param t Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string filename, IndexType t)
{
    std::string fileName = filename + "." + vectorName + "." + std::to_string(static_cast<long long>(t)) + ".mtx";
    std::cout << "snapshot for Timestep " << t << " has been written to: " << fileName << std::endl;

    vector.writeToFile(fileName);
}

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

template class KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<double>;
