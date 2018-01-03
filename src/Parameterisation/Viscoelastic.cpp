#include "Viscoelastic.hpp"
using namespace scai;

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType>::Viscoelastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(config, ctx, dist);
}

/*! \brief Initialisation with zeros
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Parameterisation::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Parameterisation::Viscoelastic<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    if (config.get<IndexType>("ModelRead")) {

        HOST_PRINT(dist->getCommunicatorPtr(), "Reading model parameter from file...\n");

        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("PartitionedIn"));
        initRelaxationMechanisms(config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency"));

        HOST_PRINT(dist->getCommunicatorPtr(), "Finished with reading of the model parameter!\n\n");

    } else {
        init(ctx, dist, config.get<ValueType>("velocityP"), config.get<ValueType>("velocityS"), config.get<ValueType>("rho"), config.get<ValueType>("tauP"), config.get<ValueType>("tauS"), config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency"));
    }

    if (config.get<IndexType>("ModelWrite")) {
        write(config.get<std::string>("ModelFilename") + ".out", config.get<IndexType>("PartitionedOut"));
        std::cout << "been here\n\n";
    }
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param velocityP_const P-wave velocity given as Scalar
 \param velocityS_const S-wave velocity given as Scalar
 \param rho_const Density given as Scalar
 \param tauP_const TauP given as Scalar
 \param tauS_const TauS given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType>::Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar velocityP_const, scai::lama::Scalar velocityS_const, scai::lama::Scalar rho_const, scai::lama::Scalar tauP_const, scai::lama::Scalar tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    init(ctx, dist, velocityP_const, velocityS_const, rho_const, tauP_const, tauS_const, numRelaxationMechanisms_in, relaxationFrequency_in);
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the five given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param sWaveModulus_const S-wave modulus given as Scalar
 \param rho_const Density given as Scalar
 \param tauP_const TauP given as Scalar
 \param tauS_const TauS given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Parameterisation::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar velocityP_const, scai::lama::Scalar velocityS_const, scai::lama::Scalar rho_const, scai::lama::Scalar tauP_const, scai::lama::Scalar tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in);
    this->initParameterisation(velocityP, ctx, dist, velocityP_const);
    this->initParameterisation(velocityS, ctx, dist, velocityS_const);
    this->initParameterisation(density, ctx, dist, rho_const);
    this->initParameterisation(tauS, ctx, dist, tauS_const);
    this->initParameterisation(tauP, ctx, dist, tauP_const);
}

/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType>::Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    init(ctx, dist, filename, partitionedIn);
}

/*! \brief Initialisator that is reading Velocity-Vector from an external files and calculates pWaveModulus
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param partitionedIn Partitioned input
 *
 *  Calculates pWaveModulus with
 */
template <typename ValueType>
void KITGPI::Parameterisation::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    std::string filenameVelocityP = filename + ".vp.mtx";
    std::string filenameVelocityS = filename + ".vs.mtx";
    std::string filenamedensity = filename + ".density.mtx";
    std::string filenameTauP = filename + ".tauP.mtx";
    std::string filenameTauS = filename + ".tauS.mtx";

    this->initParameterisation(velocityS, ctx, dist, filenameVelocityS, partitionedIn);
    this->initParameterisation(velocityP, ctx, dist, filenameVelocityP, partitionedIn);
    this->initParameterisation(density, ctx, dist, filenamedensity, partitionedIn);
    this->initParameterisation(tauS, ctx, dist, filenameTauS, partitionedIn);
    this->initParameterisation(tauP, ctx, dist, filenameTauP, partitionedIn);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType>::Viscoelastic(const Viscoelastic &rhs)
{
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
    tauS = rhs.tauS;
    tauP = rhs.tauP;
    relaxationFrequency = rhs.relaxationFrequency;
    numRelaxationMechanisms = rhs.numRelaxationMechanisms;
}

/*! \brief Write model to an external file
 *
 \param filename For the P-wave velocity ".vp.mtx" is added, for the S-wave velocity ".vs.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Parameterisation::Viscoelastic<ValueType>::write(std::string filename, IndexType partitionedOut) const
{

    std::string filenamedensity = filename + ".density.mtx";
    std::string filenameTauP = filename + ".tauP.mtx";
    std::string filenameTauS = filename + ".tauS.mtx";
    std::string filenameP = filename + ".vp.mtx";
    std::string filenameS = filename + ".vs.mtx";

    this->writeParameterisation(density, filenamedensity, partitionedOut);
    this->writeParameterisation(tauP, filenameTauP, partitionedOut);
    this->writeParameterisation(tauS, filenameTauS, partitionedOut);
    this->writeParameterisation(velocityP, filenameP, partitionedOut);
    this->writeParameterisation(velocityS, filenameS, partitionedOut);
};

/*! \brief Initialisation the relaxation mechanisms
 *
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Parameterisation::Viscoelastic<ValueType>::initRelaxationMechanisms(IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    if (numRelaxationMechanisms_in < 1) {
        COMMON_THROWEXCEPTION("The number of relaxation mechanisms should be >0 in an visco-elastic simulation")
    }
    if (relaxationFrequency_in <= 0) {
        COMMON_THROWEXCEPTION("The relaxation frequency should be >=0 in an visco-elastic simulation")
    }
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    relaxationFrequency = relaxationFrequency_in;
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType> KITGPI::Parameterisation::Viscoelastic<ValueType>::operator*(scai::lama::Scalar rhs)
{
    KITGPI::Parameterisation::Viscoelastic<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType> operator*(scai::lama::Scalar lhs, KITGPI::Parameterisation::Viscoelastic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType> &KITGPI::Parameterisation::Viscoelastic<ValueType>::operator*=(scai::lama::Scalar const &rhs)
{
    density *= rhs;
    tauS *= rhs;
    tauP *= rhs;
    velocityP *= rhs;
    velocityS *= rhs;

    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType> KITGPI::Parameterisation::Viscoelastic<ValueType>::operator+(KITGPI::Parameterisation::Viscoelastic<ValueType> const &rhs)
{
    KITGPI::Parameterisation::Viscoelastic<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType> &KITGPI::Parameterisation::Viscoelastic<ValueType>::operator+=(KITGPI::Parameterisation::Viscoelastic<ValueType> const &rhs)
{
    density += rhs.density;
    tauS += rhs.tauS;
    tauP += rhs.tauP;
    velocityP += rhs.velocityP;
    velocityS += rhs.velocityS;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType> KITGPI::Parameterisation::Viscoelastic<ValueType>::operator-(KITGPI::Parameterisation::Viscoelastic<ValueType> const &rhs)
{
    KITGPI::Parameterisation::Viscoelastic<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType> &KITGPI::Parameterisation::Viscoelastic<ValueType>::operator-=(KITGPI::Parameterisation::Viscoelastic<ValueType> const &rhs)
{
    density = density -= rhs.density;
    tauS -= rhs.tauS;
    tauP -= rhs.tauP;
    velocityP -= rhs.velocityP;
    velocityS -= rhs.velocityS;

    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Parameterisation::Viscoelastic<ValueType> &KITGPI::Parameterisation::Viscoelastic<ValueType>::operator=(KITGPI::Parameterisation::Viscoelastic<ValueType> const &rhs)
{

    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
    tauS = rhs.tauS;
    tauP = rhs.tauP;
    relaxationFrequency = rhs.relaxationFrequency;
    numRelaxationMechanisms = rhs.numRelaxationMechanisms;

    return *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is assigned.
 */
template <typename ValueType>
void KITGPI::Parameterisation::Viscoelastic<ValueType>::assign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs)
{

    density = rhs.getDensity();
    velocityP = rhs.getVelocityP();
    velocityS = rhs.getVelocityS();
    tauS = rhs.getTauS();
    tauP = rhs.getTauP();
    relaxationFrequency = rhs.getRelaxationFrequency();
    numRelaxationMechanisms = rhs.getNumRelaxationMechanisms();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is subtractet.
 */
template <typename ValueType>
void KITGPI::Parameterisation::Viscoelastic<ValueType>::minusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs)
{

    density -= rhs.getDensity();
    velocityP -= rhs.getVelocityP();
    velocityS = rhs.getVelocityS();
    tauS = rhs.getTauS();
    tauP = rhs.getTauP();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstarct model which is subtractet.
 */
template <typename ValueType>
void KITGPI::Parameterisation::Viscoelastic<ValueType>::plusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs)
{

    density += rhs.getDensity();
    velocityP += rhs.getVelocityP();
    velocityS = rhs.getVelocityS();
    tauS = rhs.getTauS();
    tauP = rhs.getTauP();
}

template class KITGPI::Parameterisation::Viscoelastic<float>;
template class KITGPI::Parameterisation::Viscoelastic<double>;
