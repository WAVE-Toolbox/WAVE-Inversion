#include "Gradient.hpp"
using namespace scai;
using namespace KITGPI;

/*! \brief Getter method for partitionedIn */
template <typename ValueType>
IndexType KITGPI::Gradient::Gradient<ValueType>::getPartitionedIn()
{
    return (PartitionedIn);
}

/*! \brief Getter method for partitionedOut */
template <typename ValueType>
IndexType KITGPI::Gradient::Gradient<ValueType>::getPartitionedOut()
{
    return (PartitionedOut);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Gradient::Gradient<ValueType>::getRelaxationFrequency() const
{
    return (relaxationFrequency);
}

/*! \brief Set relaxation frequency
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setRelaxationFrequency(ValueType const setRelaxationFrequency)
{
    relaxationFrequency = setRelaxationFrequency;
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Gradient::Gradient<ValueType>::getNumRelaxationMechanisms() const
{
    return (numRelaxationMechanisms);
}

/*! \brief Set number of Relaxation Mechanisms
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setNumRelaxationMechanisms(IndexType const setNumRelaxationMechanisms)
{
    numRelaxationMechanisms = setNumRelaxationMechanisms;
}

/*! \brief Init a single parameter by a constant value
 *
 \param vector Singel parameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param value Value which will be used to initialize the single parameter to a homogenoeus parameter vector
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::initParameterisation(scai::lama::Vector &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar value)
{

    allocateParameterisation(vector, ctx, dist);

    vector.assign(value);
}

/*! \brief Init a single parameter by reading a parameter vector from an external file
 *
 *  Reads a single parameter from an external mtx file.
 \param vector Singel parameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param filename Location of external file which will be read in
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::initParameterisation(scai::lama::Vector &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{

    allocateParameterisation(vector, ctx, dist);

    readParameterisation(vector, filename, dist, partitionedIn);

    vector.redistribute(dist);
}

/*! \brief Write singe parameter to an external file
 *
 *  Write a single parameter to an external file block.
 \param vector Single parameter which will be written to filename
 \param filename Name of file in which parameter will be written
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::writeParameterisation(scai::lama::Vector const &vector, std::string filename, IndexType partitionedOut) const
{
    PartitionedInOut::PartitionedInOut<ValueType> partitionOut;

    switch (partitionedOut) {
    case false:
        vector.writeToFile(filename);
        HOST_PRINT(vector.getDistributionPtr()->getCommunicatorPtr(), "writing " << filename << "\n");
        break;

    case true:
        partitionOut.writeToDistributedFiles(vector, filename);
        break;

    default:
        COMMON_THROWEXCEPTION("Unexpected output option!")
        break;
    }
};

/*! \brief Read a parameter from file
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::readParameterisation(scai::lama::Vector &vector, std::string filename, scai::dmemo::DistributionPtr dist, IndexType partitionedIn)
{

    PartitionedInOut::PartitionedInOut<ValueType> partitionIn;

    switch (partitionedIn) {
    case false:
        partitionIn.readFromOneFile(vector, filename, dist);
        break;

    case true:
        partitionIn.readFromDistributedFiles(vector, filename, dist);
        break;

    default:
        COMMON_THROWEXCEPTION("Unexpected input option!")
        break;
    }
};

/*! \brief Allocate a single parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::allocateParameterisation(scai::lama::Vector &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

/*! \brief Get const reference to density  parameter
 * 
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Gradient<ValueType>::getDensity()
{

    return (density);
}

/*! \brief Get const reference to density  parameter
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Gradient<ValueType>::getDensity() const
{
    return (density);
}

/*! \brief Set density  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setDensity(scai::lama::Vector const &setDensity)
{
    density = setDensity;
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Gradient<ValueType>::getVelocityP()
{
    return (velocityP);
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Gradient<ValueType>::getVelocityP() const
{
    return (velocityP);
}

/*! \brief Set P-wave velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setVelocityP(scai::lama::Vector const &setVelocityP)
{
    velocityP = setVelocityP;
}

/*! \brief Get const reference to S-wave velocity
 * 
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Gradient<ValueType>::getVelocityS()
{
    return (velocityS);
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Gradient<ValueType>::getVelocityS() const
{
    return (velocityS);
}

/*! \brief Set S-wave velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setVelocityS(scai::lama::Vector const &setVelocityS)
{
    velocityS = setVelocityS;
}

/*! \brief Get const reference to tauP
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Gradient<ValueType>::getTauP()
{
    return (tauP);
}

/*! \brief Get const reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Gradient<ValueType>::getTauP() const
{
    return (tauP);
}

/*! \brief Set tauP velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setTauP(scai::lama::Vector const &setTauP)
{
    tauP = setTauP;
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Gradient<ValueType>::getTauS()
{
    return (tauS);
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Gradient<ValueType>::getTauS() const
{
    return (tauS);
}

/*! \brief Set tauS velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setTauS(scai::lama::Vector const &setTauS)
{
    tauS = setTauS;
}

/*! \brief Overloading = Operation
 *
 \param rhs Gradient which is copied.
 */
template <typename ValueType>
KITGPI::Gradient::Gradient<ValueType> &KITGPI::Gradient::Gradient<ValueType>::operator=(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    assign(rhs);
    return *this;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Gradient<ValueType> &KITGPI::Gradient::Gradient<ValueType>::operator-=(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    minusAssign(rhs);
    return *this;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Gradient<ValueType> &KITGPI::Gradient::Gradient<ValueType>::operator+=(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    plusAssign(rhs);
    return *this;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Gradient<ValueType> &KITGPI::Gradient::Gradient<ValueType>::operator*=(scai::lama::Scalar const &rhs)
{
    timesAssign(rhs);
    return *this;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Gradient<ValueType> &KITGPI::Gradient::Gradient<ValueType>::operator*=(scai::lama::Vector const &rhs)
{
    timesAssign(rhs);
    return *this;
}

template class KITGPI::Gradient::Gradient<float>;
template class KITGPI::Gradient::Gradient<double>;
