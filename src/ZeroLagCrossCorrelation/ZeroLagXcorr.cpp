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

/*! \brief prepare for inversion
 \param setGradientType setGradientType
 \param config Configuration
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::prepareForInversion(scai::IndexType setGradientType, KITGPI::Configuration::Configuration config)
{
    gradientType = setGradientType;
    decomposeType = config.getAndCatch("decomposeType", 0);
    numRelaxationMechanisms = config.get<IndexType>("numRelaxationMechanisms");
    if (numRelaxationMechanisms > 0) {
        std::vector<ValueType> relaxationFrequency_in(numRelaxationMechanisms, 0);
        for (int l=0; l<numRelaxationMechanisms; l++) {
            if (l==0)
                relaxationFrequency_in[l] = config.get<ValueType>("relaxationFrequency");
            if (l==1)
                relaxationFrequency_in[l] = config.get<ValueType>("relaxationFrequency2");
            if (l==2)
                relaxationFrequency_in[l] = config.get<ValueType>("relaxationFrequency3");
            if (l==3)
                relaxationFrequency_in[l] = config.get<ValueType>("relaxationFrequency4");
        }
        relaxationFrequency = relaxationFrequency_in;
    }
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

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr<double>;
