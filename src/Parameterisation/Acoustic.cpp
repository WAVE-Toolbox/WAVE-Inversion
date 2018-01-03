#include "Acoustic.hpp"
#include <scai/lama/io/FileIO.hpp>
using namespace scai;
using namespace KITGPI;


/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType>::Acoustic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(config, ctx, dist);
}


/*! \brief Initialisation with zeros
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Parameterisation::Acoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
        init(ctx,dist,0.0,0.0);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Parameterisation::Acoustic<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    if (config.get<IndexType>("ModelRead")) {

        HOST_PRINT(dist->getCommunicatorPtr(), "Reading model parameter from file...\n");

         init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("PartitionedIn"));

        HOST_PRINT(dist->getCommunicatorPtr(), "Finished with reading of the model parameter!\n\n");

    } else {
        init(ctx, dist, config.get<ValueType>("velocityP"), config.get<ValueType>("rho"));
    }

    if (config.get<IndexType>("ModelWrite")) {
        write(config.get<std::string>("ModelFilename") + ".out", config.get<IndexType>("PartitionedOut"));
    }
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param rho_const Density given as Scalar
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType>::Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus_const, scai::lama::Scalar rho_const)
{
    init(ctx, dist, pWaveModulus_const, rho_const);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param rho_const Density given as Scalar
 */
template <typename ValueType>
void KITGPI::Parameterisation::Acoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar velocityP_const, scai::lama::Scalar rho_const)
{
    this->initParameterisation(velocityP, ctx, dist, velocityP_const);
    this->initParameterisation(density, ctx, dist, rho_const);
}

/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added and for density ".density.mtx" is added.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType>::Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    init(ctx, dist, filename, partitionedIn);
}

/*! \brief Initialisator that is reading Velocity-Vector
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the first Velocity-Vector "filename".vp.mtx" is added and for density "filename+".density.mtx" is added.
 \param partitionedIn Partitioned input
 *
 */
template <typename ValueType>
void KITGPI::Parameterisation::Acoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    std::string filenameVelocityP = filename + ".vp.mtx";
    std::string filenamedensity = filename + ".density.mtx";

    this->initParameterisation(velocityP, ctx, dist, filenameVelocityP, partitionedIn);
    this->initParameterisation(density, ctx, dist, filenamedensity, partitionedIn);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType>::Acoustic(const Acoustic &rhs)
{
    velocityP = rhs.velocityP;
    density = rhs.density;
}

/*! \brief Write model to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added and for density ".density.mtx" is added.
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Parameterisation::Acoustic<ValueType>::write(std::string filename, IndexType partitionedOut) const
{
    std::string filenameP = filename + ".vp.mtx";
    std::string filenamedensity = filename + ".density.mtx";
    
    this->writeParameterisation(density, filenamedensity, partitionedOut);
    this->writeParameterisation(velocityP, filenameP, partitionedOut);
};





/*! \brief Get reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Parameterisation::Acoustic<ValueType>::getVelocityS()
{
    COMMON_THROWEXCEPTION("The S-wave velocity is not defined in an acoustic simulation.")
    return (velocityS);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Parameterisation::Acoustic<ValueType>::getTauP()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Parameterisation::Acoustic<ValueType>::getTauS()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauS);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Parameterisation::Acoustic<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an elastic modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Parameterisation::Acoustic<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an elastic modelling")
    return (numRelaxationMechanisms);
}



/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType> KITGPI::Parameterisation::Acoustic<ValueType>::operator*(scai::lama::Scalar rhs)
{
    KITGPI::Parameterisation::Acoustic<ValueType> result(*this);
    result *= rhs;
    return result;   
}

/*! \brief non-member function to multiply (scalar as left operand)
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType> operator*(scai::lama::Scalar lhs, KITGPI::Parameterisation::Acoustic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType> &KITGPI::Parameterisation::Acoustic<ValueType>::operator*=(scai::lama::Scalar const &rhs)
{
    density *= rhs;
    velocityP *= rhs;
    
    return *this;

}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType> KITGPI::Parameterisation::Acoustic<ValueType>::operator+(KITGPI::Parameterisation::Acoustic<ValueType> const &rhs)
{
    KITGPI::Parameterisation::Acoustic<ValueType> result(*this);
    result += rhs;
    return result;    
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType> &KITGPI::Parameterisation::Acoustic<ValueType>::operator+=(KITGPI::Parameterisation::Acoustic<ValueType> const &rhs)
{
    density += rhs.density;
    velocityP += rhs.velocityP;

        return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType> KITGPI::Parameterisation::Acoustic<ValueType>::operator-(KITGPI::Parameterisation::Acoustic<ValueType> const &rhs)
{
    KITGPI::Parameterisation::Acoustic<ValueType> result(*this);
    result -= rhs;
    return result; 
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType> &KITGPI::Parameterisation::Acoustic<ValueType>::operator-=(KITGPI::Parameterisation::Acoustic<ValueType> const &rhs)
{
     density -= rhs.density;
    velocityP -= rhs.velocityP;
        return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Parameterisation::Acoustic<ValueType> &KITGPI::Parameterisation::Acoustic<ValueType>::operator=(KITGPI::Parameterisation::Acoustic<ValueType> const &rhs)
{
    // why does rhs.density not work (density = protected) 
    velocityP = rhs.velocityP;
    density = rhs.density;
    return *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is assigned.
 */
template <typename ValueType>
void KITGPI::Parameterisation::Acoustic<ValueType>::assign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs)
{

         density = rhs.getDensity();
        velocityP = rhs.getVelocityP();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is subtractet.
 */
template <typename ValueType>
void KITGPI::Parameterisation::Acoustic<ValueType>::minusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs)
{

         density -= rhs.getDensity();
        velocityP -= rhs.getVelocityP();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstarct model which is subtractet.
 */
template <typename ValueType>
void KITGPI::Parameterisation::Acoustic<ValueType>::plusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs)
{

         density += rhs.getDensity();
        velocityP += rhs.getVelocityP();
}


template class KITGPI::Parameterisation::Acoustic<double>;
template class KITGPI::Parameterisation::Acoustic<float>;
