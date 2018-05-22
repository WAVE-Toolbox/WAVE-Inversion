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

    vector = lama::zero<lama::DenseVector<ValueType>>( dist, ctx );

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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string type, IndexType t)
{
    std::string fileName = "wavefields/wavefield" + type + "." + vectorName + "." + std::to_string(static_cast<long long>(t)) + ".mtx";
    std::cout << "snapshot for Timestep " << t << "has been written to: " << fileName;

    vector.writeToFile(fileName);
}

//! \brief Getter routine for VSum wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::getVSum() const
{
    return (VSum);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::getP() const
{
    return (P);
}

//! \brief Getter routine for ShearStress
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::getShearStress() const
{
    return (ShearStress);
}

//! \brief Getter routine for NormalStressDiff
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::getNormalStressDiff() const
{
    return (NormalStressDiff);
}

//! \brief Getter routine for NormalStressSum
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::getNormalStressSum() const
{
    return (NormalStressSum);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr<double>;
