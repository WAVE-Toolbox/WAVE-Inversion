#include "Elastic.hpp"

using namespace scai;

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType>::Elastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(config, ctx, dist,0.0, 0.0, 0.0);
}

/*! \brief Initialisation with zeros
 *
  \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::init(Configuration::Configuration const &config,scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(config,ctx, dist, 0.0, 0.0, 0.0);
}


/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param sWaveModulus_const S-wave modulus given as Scalar
 \param rho Density given as Scalar
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::init(Configuration::Configuration const &config,scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar velocityP_const, scai::lama::Scalar velocityS_const, scai::lama::Scalar rho)
{
    invertForVp=config.get<bool>("invertForVp");
    invertForVs=config.get<bool>("invertForVs");
    invertForDensity=config.get<bool>("invertForDensity");
    this->initParameterisation(velocityP, ctx, dist, velocityP_const);
    this->initParameterisation(velocityS, ctx, dist, velocityS_const);
    this->initParameterisation(density, ctx, dist, rho);
}



//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType>::Elastic(const Elastic &rhs)
{
    invertForVp=rhs.invertForVp;
    invertForVs=rhs.invertForVs;
    invertForDensity=rhs.invertForDensity;
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
}

/*! \brief Write model to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::write(std::string filename, IndexType partitionedOut) const
{
    std::string filenameP = filename + ".vp.mtx";
    std::string filenameS = filename + ".vs.mtx";
    std::string filenamedensity = filename + ".density.mtx";

    this->writeParameterisation(density, filenamedensity, partitionedOut);
    this->writeParameterisation(velocityP, filenameP, partitionedOut);
    this->writeParameterisation(velocityS, filenameS, partitionedOut);
};

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Elastic<ValueType>::getTauP()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Gradient::Elastic<ValueType>::getTauS()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauS);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Gradient::Elastic<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an elastic modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Gradient::Elastic<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an elastic modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> KITGPI::Gradient::Elastic<ValueType>::operator*(scai::lama::Scalar rhs)
{
    KITGPI::Gradient::Elastic<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> operator*(scai::lama::Scalar lhs, KITGPI::Gradient::Elastic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> &KITGPI::Gradient::Elastic<ValueType>::operator*=(scai::lama::Scalar const &rhs)
{
    density *= rhs;
    velocityP *= rhs;
    velocityS *= rhs;

    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> KITGPI::Gradient::Elastic<ValueType>::operator+(KITGPI::Gradient::Elastic<ValueType> const &rhs)
{
    KITGPI::Gradient::Elastic<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> &KITGPI::Gradient::Elastic<ValueType>::operator+=(KITGPI::Gradient::Elastic<ValueType> const &rhs)
{
    density += rhs.density;
    velocityP += rhs.velocityP;
    velocityS += rhs.velocityS;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> KITGPI::Gradient::Elastic<ValueType>::operator-(KITGPI::Gradient::Elastic<ValueType> const &rhs)
{
    KITGPI::Gradient::Elastic<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> &KITGPI::Gradient::Elastic<ValueType>::operator-=(KITGPI::Gradient::Elastic<ValueType> const &rhs)
{
    density -= rhs.density;
    velocityP -= rhs.velocityP;
    velocityS -= rhs.velocityS;

    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Gradient which is copied.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> &KITGPI::Gradient::Elastic<ValueType>::operator=(KITGPI::Gradient::Elastic<ValueType> const &rhs)
{
    invertForVp=rhs.invertForVp;
    invertForVs=rhs.invertForVs;
    invertForDensity=rhs.invertForDensity;	
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;

    return *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::assign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    invertForVp=rhs.invertForVp;
    invertForVs=rhs.invertForVs;
    invertForDensity=rhs.invertForDensity;
    density = rhs.getDensity();
    velocityP = rhs.getVelocityP();
    velocityS = rhs.getVelocityS();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{

    density -= rhs.getDensity();
    velocityP -= rhs.getVelocityP();
    velocityS -= rhs.getVelocityS();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{

    density += rhs.getDensity();
    velocityP += rhs.getVelocityP();
    velocityS += rhs.getVelocityS();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::timesAssign(scai::lama::Scalar const &rhs)
{
    density *= rhs;
    velocityP *= rhs;
    velocityS *= rhs;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract model.
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs){
    scai::lama::DenseVector<ValueType> temp;
    temp = lhs.getVelocityP() - rhs.getVelocityP();
    lhs.setVelocityP(temp);
    temp = lhs.getVelocityS() - rhs.getVelocityS();
    lhs.setVelocityS(temp);
    temp = lhs.getDensity() - rhs.getDensity();
    lhs.setDensity(temp);
    
};

/*! \brief function for scaling the gradients with the model parameter 
 *
 \param model Abstract model.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    if (invertForVp) {
        velocityP *= 1 / velocityP.maxNorm() * model.getVelocityP().maxNorm();
    }
    
    if (invertForVs) {
        velocityS *= 1 / velocityS.maxNorm() * model.getVelocityS().maxNorm();
    }   
    
    if (invertForDensity) {
        density *= 1 / density.maxNorm() * model.getDensity().maxNorm();
    }
}

template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT)
{
    // Only implementedfor 2D!
    //dt should be in cross correlation!

   // Lambda and Mu gradients
	
    scai::lama::DenseVector<ValueType> gradLambda;
    scai::lama::DenseVector<ValueType> gradMu;  
    scai::lama::DenseVector<ValueType> temp;
    
    
    gradLambda = model.getVelocityP();
    gradLambda *= model.getVelocityP();
    temp = model.getVelocityS();
    temp *= model.getVelocityS();
    gradLambda-=temp;
    gradLambda*=gradLambda;
    gradLambda*=model.getDensity();
    gradLambda*=model.getDensity();
    gradLambda*=4;
    gradLambda.invert();
    
    gradLambda *= correlatedWavefields.getNormalStressSum();
    gradLambda *= -DT;

    if ((invertForVs) || (invertForDensity)){
	  gradMu=gradLambda;
 
	  temp = model.getVelocityS();
          temp *= model.getVelocityS();
	  temp *= model.getDensity();
	  temp *= temp;
	  temp*=4;
	  temp.invert();
	  temp *= correlatedWavefields.getNormalStressDiff();
	  temp *= -DT;
	  
	  gradMu += temp;
	
	  temp = model.getVelocityS();
          temp *= model.getVelocityS();
	  temp *= model.getDensity();
	  temp *= temp;
	  temp.invert();
	  temp *= correlatedWavefields.getShearStress();
	  temp *= -DT;
	  
	  gradMu += temp;
    }
    
    // vp, vs , rho gradients
    
    if (invertForVp) {
	//grad_vp = 2*rho*vp*grad_bulk
        velocityP = 2 * gradLambda;
        velocityP *= model.getDensity();
        velocityP *= model.getVelocityP();
    }

    if (invertForVs) {

        velocityS = -4 * gradLambda;
        velocityS *= model.getDensity();
        velocityS *= model.getVelocityS();
	
	temp=2*gradMu;
	temp*=model.getDensity();
	temp*=model.getVelocityS();
	
	velocityS+=temp;
    }
    
    if (invertForDensity) {

        density = model.getVelocityP();
        density *= model.getVelocityP();
        temp = 2 * model.getVelocityS();
        temp *= model.getVelocityS();
        density-=temp;
	density*=gradLambda;
	
	temp = model.getVelocityS();
	temp*= model.getVelocityS();
	temp*= gradMu;
	
	density+=temp;
	
        density += DT * correlatedWavefields.getVSum();
    }
}

template class KITGPI::Gradient::Elastic<float>;
template class KITGPI::Gradient::Elastic<double>;
