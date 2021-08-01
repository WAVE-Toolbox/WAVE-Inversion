#include "ZeroLagXcorr.hpp"

using namespace scai;

/*! \brief Reset a single wavefield to zero.
 \param vector Vector to be reset to 0
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::resetWavefield(scai::lama::DenseVector<ValueType> &vector)
{
    vector = 0;
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
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
 \param type Wavefield-type (acoustic, elastic, viscoelastic)
 \param t Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string filename, IndexType t)
{
    std::string fileName = filename + "." + vectorName + "." + std::to_string(static_cast<long long>(t)) + ".mtx";
    std::cout << "snapshot for Timestep " << t << " has been written to: " << fileName << std::endl;

    vector.writeToFile(fileName);
}

//! \brief Getter routine for VSum wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::getXcorrRho() const
{
    return (xcorrRho);
}

//! \brief Getter routine for XcorrLambda
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::getXcorrLambda() const
{
    return (xcorrLambda);
}

//! \brief Getter routine for xcorrMuA
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::getXcorrMuA() const
{
    return (xcorrMuA);
}

//! \brief Getter routine for xcorrMuB
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::getXcorrMuB() const
{
    return (xcorrMuB);
}

//! \brief Getter routine for xcorrMuC
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::getXcorrMuC() const
{
    return (xcorrMuC);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr<double>;
