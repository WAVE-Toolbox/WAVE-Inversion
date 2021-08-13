#include "Viscoelastic.hpp"
#include <scai/dmemo/SingleDistribution.hpp>

using namespace scai;

/*! \brief Constructor that is using the Configuration class
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Gradient::Viscoelastic<ValueType>::Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "viscoelastic";
    init(ctx, dist, 0.0, 0.0, 0.0);
}

/*! \brief Initialisation with zeros
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist, 0.0, 0.0, 0.0);
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
void KITGPI::Gradient::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho_const)
{
    this->initParameterisation(velocityP, ctx, dist, velocityP_const);
    this->initParameterisation(velocityS, ctx, dist, velocityS_const);
    this->initParameterisation(density, ctx, dist, rho_const);
    this->initParameterisation(porosity, ctx, dist, 0.0);
    this->initParameterisation(saturation, ctx, dist, 0.0);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::Viscoelastic<ValueType>::Viscoelastic(const Viscoelastic &rhs)
{
    equationType = rhs.equationType;
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
}

/*! \brief Set all parameter to zero.
*/
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::resetGradient()
{
    this->resetParameter(velocityP);
    this->resetParameter(velocityS);
    this->resetParameter(density);
    this->resetParameter(porosity);
    this->resetParameter(saturation);
}
/*! \brief Write model to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 \param fileFormat format of output file
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::write(std::string filename, IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const
{
    if (workflow.getInvertForVp()) {
        std::string filenameP = filename + ".vp";
        this->writeParameterisation(velocityP, filenameP, fileFormat);
    }

    if (workflow.getInvertForVs()) {
        std::string filenameS = filename + ".vs";
        this->writeParameterisation(velocityS, filenameS, fileFormat);
    }

    if (workflow.getInvertForDensity()) {
        std::string filenamedensity = filename + ".density";
        this->writeParameterisation(density, filenamedensity, fileFormat);
    }
    
    if (workflow.getInvertForPorosity()) { 
        std::string filenamePorosity = filename + ".porosity";
        this->writeParameterisation(porosity, filenamePorosity, fileFormat);
    }
        
    if (workflow.getInvertForSaturation()) { 
        std::string filenameSaturation = filename + ".saturation";
        this->writeParameterisation(saturation, filenameSaturation, fileFormat);
    }
};

/*! \brief Get equationType (viscoelastic)
 */
template <typename ValueType>
std::string KITGPI::Gradient::Viscoelastic<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Viscoelastic<ValueType>::getTauP()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an viscoelastic modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Viscoelastic<ValueType>::getTauS()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an viscoelastic modelling")
    return (tauS);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Gradient::Viscoelastic<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an viscoelastic modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Gradient::Viscoelastic<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an viscoelastic modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::Viscoelastic<ValueType> KITGPI::Gradient::Viscoelastic<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Gradient::Viscoelastic<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief Free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Gradient::Viscoelastic<ValueType> operator*(ValueType lhs, KITGPI::Gradient::Viscoelastic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::Viscoelastic<ValueType> &KITGPI::Gradient::Viscoelastic<ValueType>::operator*=(ValueType const &rhs)
{
    if (workflowInner.getInvertForDensity()) 
        density *= rhs;
    if (workflowInner.getInvertForVp()) 
        velocityP *= rhs;
    if (workflowInner.getInvertForVs()) 
        velocityS *= rhs;
    if (workflowInner.getInvertForPorosity()) 
        porosity *= rhs;
    if (workflowInner.getInvertForSaturation()) 
        saturation *= rhs;

    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::Viscoelastic<ValueType> KITGPI::Gradient::Viscoelastic<ValueType>::operator+(KITGPI::Gradient::Viscoelastic<ValueType> const &rhs)
{
    KITGPI::Gradient::Viscoelastic<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::Viscoelastic<ValueType> &KITGPI::Gradient::Viscoelastic<ValueType>::operator+=(KITGPI::Gradient::Viscoelastic<ValueType> const &rhs)
{
    density += rhs.density;
    velocityP += rhs.velocityP;
    velocityS += rhs.velocityS;
    porosity += rhs.porosity;
    saturation += rhs.saturation;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Viscoelastic<ValueType> KITGPI::Gradient::Viscoelastic<ValueType>::operator-(KITGPI::Gradient::Viscoelastic<ValueType> const &rhs)
{
    KITGPI::Gradient::Viscoelastic<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Viscoelastic<ValueType> &KITGPI::Gradient::Viscoelastic<ValueType>::operator-=(KITGPI::Gradient::Viscoelastic<ValueType> const &rhs)
{
    density -= rhs.density;
    velocityP -= rhs.velocityP;
    velocityS -= rhs.velocityS;
    porosity -= rhs.porosity;
    saturation -= rhs.saturation;

    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Gradient which is copied.
 */
template <typename ValueType>
KITGPI::Gradient::Viscoelastic<ValueType> &KITGPI::Gradient::Viscoelastic<ValueType>::operator=(KITGPI::Gradient::Viscoelastic<ValueType> const &rhs)
{
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
    porosity = rhs.porosity;
    saturation = rhs.saturation;

    return *this;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::assign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    density = rhs.getDensity();
    velocityP = rhs.getVelocityP();
    velocityS = rhs.getVelocityS();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{

    density -= rhs.getDensity();
    velocityP -= rhs.getVelocityP();
    velocityS -= rhs.getVelocityS();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
}

/*! \brief Function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{

    density += rhs.getDensity();
    velocityP += rhs.getVelocityP();
    velocityS += rhs.getVelocityS();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::timesAssign(ValueType const &rhs)
{
    if (workflowInner.getInvertForDensity()) 
        density *= rhs;
    if (workflowInner.getInvertForVp()) 
        velocityP *= rhs;
    if (workflowInner.getInvertForVs()) 
        velocityS *= rhs;
    if (workflowInner.getInvertForPorosity()) 
        porosity *= rhs;
    if (workflowInner.getInvertForSaturation()) 
        saturation *= rhs;
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::timesAssign(scai::lama::Vector<ValueType> const &rhs)
{
    density *= rhs;
    velocityP *= rhs;
    velocityS *= rhs;
    porosity *= rhs;
    saturation *= rhs;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract model.
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs)
{    
    scai::lama::DenseVector<ValueType> temp;   
    if (lhs.getParameterisation() == 1 || lhs.getParameterisation() == 2) {   
        if (workflowInner.getInvertForPorosity()) {
            temp = lhs.getPorosity() - rhs.getPorosity();
            lhs.setPorosity(temp);
        }
        if (workflowInner.getInvertForSaturation()) {
            temp = lhs.getSaturation() - rhs.getSaturation();
            lhs.setSaturation(temp);     
        }
    } else if (lhs.getParameterisation() == 3) {  
        if (workflowInner.getInvertForVp()) {
            temp = lhs.getVelocityP() - rhs.getVelocityP();
            lhs.setVelocityP(temp);
        }  
        if (workflowInner.getInvertForVs()) {
            temp = lhs.getVelocityS() - rhs.getVelocityS();
            lhs.setVelocityS(temp);
        }  
        if (workflowInner.getInvertForDensity()) {
            temp = lhs.getDensity() - rhs.getDensity();
            lhs.setDensity(temp);
        }
    } else if (lhs.getParameterisation() == 0) {
        scai::lama::DenseVector<ValueType> rho;
        scai::lama::DenseVector<ValueType> mu;
        rho = lhs.getDensity();
        mu = scai::lama::pow(lhs.getVelocityS(), 2);
        mu *= rho;
        if (workflowInner.getInvertForDensity()) {
            temp = lhs.getDensity() - rhs.getDensity();
            lhs.setDensity(temp);
        }        
        if (workflowInner.getInvertForVs()) {
            temp = mu - rhs.getVelocityS();
            temp /= rho;
            temp = scai::lama::sqrt(temp);
            Common::replaceInvalid<ValueType>(temp, 0.0);
            lhs.setVelocityS(temp);
        }
        if (workflowInner.getInvertForVp()) {
            scai::lama::DenseVector<ValueType> lambda;
            lambda = scai::lama::pow(lhs.getVelocityP(), 2);
            lambda *= rho;
            lambda -= 2 * mu;                    
            lambda -= rhs.getVelocityP();             
            lambda += 2 * mu;
            lambda /= rho;
            lambda = scai::lama::sqrt(lambda);
            Common::replaceInvalid<ValueType>(lambda, 0.0);
            lhs.setVelocityP(lambda);
        }
    }
};

/*! \brief Function for summing the gradients of all shot domains
 *
 \param commInterShot inter shot communication pointer
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot)
{
    /*reduction between shot domains.
    each shot domain may have a different distribution of (gradient) vectors. This happens if geographer is used (different result for dist on each shot domain even for homogenous architecture) or on heterogenous architecture. In this case even the number of processes on each domain can vary. Therfore it is necessary that only one process per shot domain communicates all data.
    */
    
    //get information from distributed vector
    auto size = velocityP.size();
    auto dist = velocityP.getDistributionPtr();
    auto comm = dist->getCommunicatorPtr();
    
    // create single distribution, only master process owns the complete vector (no distribution).
    
    int shotMaster = 0;
    auto singleDist = std::make_shared<dmemo::SingleDistribution>( size, comm, shotMaster );
    
    //redistribute vector to master process
    // (this may cause memory issues for big models)
    velocityP.redistribute(singleDist);
    
    //reduce local array (size of local array is !=0 only for master process)
    commInterShot->sumArray(velocityP.getLocalValues());
    
    //redistribute vector to former partition
    velocityP.redistribute(dist);
    
    velocityS.redistribute(singleDist);
    commInterShot->sumArray(velocityS.getLocalValues());
    velocityS.redistribute(dist);
    
    density.redistribute(singleDist);
    commInterShot->sumArray(density.getLocalValues());
    density.redistribute(dist);
    
    porosity.redistribute(singleDist);
    commInterShot->sumArray(porosity.getLocalValues());
    porosity.redistribute(dist);
    
    saturation.redistribute(singleDist);
    commInterShot->sumArray(saturation.getLocalValues());
    saturation.redistribute(dist);
}

/*! \brief If stream configuration is used, set a gradient into the big gradient
 \param model model
 \param gradientPerShot gradient per shot
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::sumGradientPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth)
{
    auto distBig = density.getDistributionPtr();
    auto dist = gradientPerShot.getDensity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> recoverMatrix = model.getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd));
    recoverMatrix.assignTranspose(recoverMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = model.getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd), boundaryWidth);
    scai::lama::SparseVector<ValueType> recoverVector;
    recoverVector = 1.0 - eraseVector;
    
    scai::lama::DenseVector<ValueType> temp;
    
    temp = recoverMatrix * gradientPerShot.getVelocityP(); //transform pershot into big model
    temp *= recoverVector;
    velocityP *= eraseVector;
    velocityP += temp; //take over the values
  
    temp = recoverMatrix * gradientPerShot.getVelocityS(); //transform pershot into big model
    temp *= recoverVector;
    velocityS *= eraseVector;
    velocityS += temp; //take over the values

    temp = recoverMatrix * gradientPerShot.getTauP(); //transform pershot into big model
    temp *= recoverVector;
    tauP *= eraseVector;
    tauP += temp; //take over the values
  
    temp = recoverMatrix * gradientPerShot.getTauS(); //transform pershot into big model
    temp *= recoverVector;
    tauS *= eraseVector;
    tauS += temp; //take over the values
    
    temp = recoverMatrix * gradientPerShot.getDensity(); //transform pershot into big model
    temp *= recoverVector;
    density *= eraseVector;
    density += temp; //take over the values
    
    temp = recoverMatrix * gradientPerShot.getPorosity(); //transform pershot into big model
    temp *= recoverVector;
    porosity *= eraseVector;
    porosity += temp; //take over the values
    
    temp = recoverMatrix * gradientPerShot.getSaturation(); //transform pershot into big model
    temp *= recoverVector;
    saturation *= eraseVector;
    saturation += temp; //take over the values
}

/*! \brief Function for scaling the gradients with the model parameter 
 *
 \param model Abstract model.
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{
    ValueType maxValue = 1;      
    
    IndexType scaleGradient = config.get<IndexType>("scaleGradient");
    if (workflow.getInvertForVp() && velocityP.maxNorm() != 0) {
        if (scaleGradient == 1) {
            if (model.getParameterisation() == 3) {
                maxValue = model.getVelocityP().maxNorm();
            } else if (model.getParameterisation() == 0) {
                scai::lama::DenseVector<ValueType> lambda;
                scai::lama::DenseVector<ValueType> mu;
                lambda = scai::lama::pow(model.getVelocityP(), 2);
                lambda *= model.getDensity();
                mu = scai::lama::pow(model.getVelocityS(), 2);
                mu *= model.getDensity();
                lambda -= 2 * mu;
                maxValue = lambda.maxNorm();
            }
        } else if (scaleGradient == 2) {
            if (model.getParameterisation() == 3) {
                maxValue = config.get<ValueType>("upperVPTh") - config.get<ValueType>("lowerVPTh");
            } else if (model.getParameterisation() == 0) {
                ValueType lambdaMax;
                ValueType lambdaMin;
                ValueType muMax;
                ValueType muMin;
                lambdaMax = pow(config.get<ValueType>("upperVPTh"), 2);
                lambdaMax *= config.get<ValueType>("upperDensityTh");
                lambdaMin = pow(config.get<ValueType>("lowerVPTh"), 2);
                lambdaMin *= config.get<ValueType>("lowerDensityTh");
                muMax = pow(config.get<ValueType>("upperVSTh"), 2);
                muMax *= config.get<ValueType>("upperDensityTh");
                muMin = pow(config.get<ValueType>("lowerVSTh"), 2);
                muMin *= config.get<ValueType>("lowerDensityTh");
                lambdaMax -= 2 * muMin;
                lambdaMin -= 2 * muMax;
                maxValue = lambdaMax - lambdaMin;
            }
        }
        velocityP *= 1 / velocityP.maxNorm() * maxValue;
    }
    
    if (workflow.getInvertForVs() && velocityS.maxNorm() != 0) {
        if (scaleGradient == 1) {
            if (model.getParameterisation() == 3) {
                maxValue = model.getVelocityS().maxNorm();
            } else if (model.getParameterisation() == 0) {
                scai::lama::DenseVector<ValueType> mu;
                mu = scai::lama::pow(model.getVelocityS(), 2);
                mu *= model.getDensity();
                maxValue = mu.maxNorm();
            }
        } else if (scaleGradient == 2) {
            if (model.getParameterisation() == 3) {
                maxValue = config.get<ValueType>("upperVSTh") - config.get<ValueType>("lowerVSTh");
            } else if (model.getParameterisation() == 0) {
                ValueType muMax;
                ValueType muMin;
                muMax = pow(config.get<ValueType>("upperVSTh"), 2);
                muMax *= config.get<ValueType>("upperDensityTh");
                muMin = pow(config.get<ValueType>("lowerVSTh"), 2);
                muMin *= config.get<ValueType>("lowerDensityTh");
                maxValue = muMax - muMin;
            }
        }
        velocityS *= 1 / velocityS.maxNorm() * maxValue;
    }
    
    if (workflow.getInvertForDensity() && density.maxNorm() != 0) {
        if (model.getParameterisation() != 1 && model.getParameterisation() != 2) {
            if (scaleGradient == 1) {
                maxValue = model.getDensity().maxNorm();
            } else if (scaleGradient == 2) {
                maxValue = config.get<ValueType>("upperDensityTh") - config.get<ValueType>("lowerDensityTh");
            }
        }
        density *= 1 / density.maxNorm() * maxValue;
    }    
    
    if (workflow.getInvertForPorosity() && porosity.maxNorm() != 0) {
        if (model.getParameterisation() == 1 || model.getParameterisation() == 2) {
            if (scaleGradient == 1) {
                maxValue = model.getPorosity().maxNorm();
            } else if (scaleGradient == 2) {
                maxValue = config.getAndCatch("upperPorosityTh", 1.0) - config.getAndCatch("lowerPorosityTh", 0.0);
            }
        }        
        porosity *= 1 / porosity.maxNorm() * maxValue;
    }    
    
    if (workflow.getInvertForSaturation() && saturation.maxNorm() != 0) { 
        if (model.getParameterisation() == 1 || model.getParameterisation() == 2) {
            if (scaleGradient == 1) {
                maxValue = model.getSaturation().maxNorm();
            } else if (scaleGradient == 2) {
                maxValue = config.getAndCatch("upperSaturationTh", 1.0) - config.getAndCatch("lowerSaturationTh", 0.0);
            }
        }              
        saturation *= 1 / saturation.maxNorm() * maxValue;        
    }
}

/*! \brief Function for normalizing the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::normalize()
{
    if (this->getNormalizeGradient()) {
        ValueType gradientMax = velocityP.maxNorm();
        if (gradientMax != 0)
            velocityP *= 1 / gradientMax;
        gradientMax = velocityS.maxNorm();
        if (gradientMax != 0)
            velocityS *= 1 / gradientMax;
        gradientMax = density.maxNorm();
        if (gradientMax != 0)
            density *= 1 / gradientMax;
        gradientMax = porosity.maxNorm();
        if (gradientMax != 0)
            porosity *= 1 / gradientMax;
        gradientMax = saturation.maxNorm();
        if (gradientMax != 0)
            saturation *= 1 / gradientMax;
    }
}

/*! \brief Function for calculating the viscoelastic gradients from the cross correlation and the model parameter 
 *
 \param model Abstract model.
 \param correlatedWavefields Abstract xCorr.
 \param DT Temporal discretization 
 \param workflow 
 *
 \f{eqnarray*}
   \nabla_{\lambda} E &=& - \mathrm{d}t \frac{1}{(N \lambda+2\mu)^2} \cdot X_{\lambda} \\
   \nabla_{\mu} E &=& \mathrm{d}t \left[ \frac{N \lambda^2 + 4 \mu \lambda}{2 \mu^2 (N \lambda+2\mu)^2} - \frac{1}{2 \mu^2} \right] \cdot X_{\mu,A} + \mathrm{d}t \frac{N \lambda^2 + 4 \mu \lambda}{2 \mu^2 (N \lambda+2\mu)^2} \cdot X_{\mu,B} - \mathrm{d}t  \frac{1}{\mu^2} \cdot X_{\mu,C} \\ \\
   \nabla_{v_{\mathrm{p}}} E &=& 2 \rho v_{\mathrm{p}} \cdot \nabla_{\lambda} E \\
   \nabla_{v_{\mathrm{s}}} E &=& 2 \rho v_{\mathrm{s}} \cdot \nabla_{\mu} E - 4 \rho v_{\mathrm{s}} \nabla_{\lambda} E \\
   \nabla_{\rho} E &=& ( v_{\mathrm{p}}^2-2v_{\mathrm{s}}^2 ) \cdot \nabla_{\lambda} E + v_{\mathrm{s}}^2 \cdot \nabla_{\mu} E - \mathrm{d}t \cdot X_{\rho}
 \f}
 *
 * with \f$ N \f$ as the number of dimensions, \f$ \mu = \rho v_{\mathrm{s}}^2 \f$ and \f$ \lambda = \rho v_{\mathrm{p}}^2 - 2\mu\f$.
 *
 \sa{KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::update} for the cross-correlations \f$ (X_{\lambda},X_{\rho},X_{\mu,A},X_{\mu,B},X_{\mu,C}) \f$ in 2D 
 \sa{KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<ValueType>::update} for the cross-correlations \f$ (X_{\lambda},X_{\rho},X_{\mu,A},X_{\mu,B},X_{\mu,C}) \f$ in 3D
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //dt should be in cross correlation!

    scai::lama::DenseVector<ValueType> gradLambda;
    scai::lama::DenseVector<ValueType> gradMu;
    scai::lama::DenseVector<ValueType> gradRho;
    scai::lama::DenseVector<ValueType> gradK_sat;
    scai::lama::DenseVector<ValueType> temp;
    scai::lama::DenseVector<ValueType> temp2;
    scai::lama::DenseVector<ValueType> lambda;
    scai::lama::DenseVector<ValueType> mu;

    //dimension
    int N = correlatedWavefields.getNumDimension();

    mu = scai::lama::pow(model.getVelocityS(), 2);
    mu *= model.getDensity();

    lambda = scai::lama::pow(model.getVelocityP(), 2);
    lambda *= model.getDensity();
    lambda -= 2 * mu;

    // Lambda and Mu gradients
    //-1/(N*lambda+2mu)^2)
    gradLambda = N * lambda;
    gradLambda += 2 * mu;
    gradLambda = scai::lama::pow(gradLambda, -2);

    gradLambda *= correlatedWavefields.getXcorrLambda();
    gradLambda *= -DT;
    Common::replaceInvalid<ValueType>(gradLambda, 0.0);

    scai::hmemo::ContextPtr ctx = gradLambda.getContextPtr();
    scai::dmemo::DistributionPtr dist = gradLambda.getDistributionPtr();

    if (workflow.getInvertForVs() || workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        //(N*lambda^2+4mu*lambda)/(2mu^2(N*lambda+2mu)^2)

        //temp2=>B
        temp2 = scai::lama::pow(lambda, 2);
        temp2 *= N;
        temp = 4 * mu;
        temp *= lambda;
        temp2 += temp;

        temp = scai::lama::pow(mu, 2);
        temp *= 2;
        temp2 /= temp;

        temp = N * lambda;
        temp += 2 * mu;
        temp = scai::lama::pow(temp, 2);
        temp2 /= temp;

        gradMu = temp2;
        gradMu *= correlatedWavefields.getXcorrMuB();

        //temp2=>A
        temp = scai::lama::pow(mu, -2);
        temp /= -2;
        temp2 += temp;

        temp2 *= correlatedWavefields.getXcorrMuA();
        gradMu += temp2;

        //temp2=>C
        temp2 = mu;
        temp2 *= mu;
        temp2 = 1 / temp2;
        temp2 *= correlatedWavefields.getXcorrMuC();
        gradMu -= temp2;

        gradMu *= DT;
        Common::replaceInvalid<ValueType>(gradMu, 0.0);
                
        gradK_sat = gradLambda;
        gradK_sat += 3 / 2 * gradMu; 
        
        gradRho = -DT * correlatedWavefields.getXcorrRho();         
    }

    // vp, vs , rho gradients
    if (workflow.getInvertForVp()) {
        if (model.getParameterisation() == 0) {
            velocityP = gradLambda;
        } else if (model.getParameterisation() == 3) {
            //grad_vp = 2*rho*vp*grad_lambda
            velocityP = 2 * gradLambda;
            velocityP *= model.getDensity();
            velocityP *= model.getVelocityP();
        }
    } else {
        this->initParameterisation(velocityP, ctx, dist, 0.0);
    }

    if (workflow.getInvertForVs()) {
        if (model.getParameterisation() == 0) {
            velocityS = gradMu;
        } else if (model.getParameterisation() == 3) {
            velocityS = -4 * gradLambda;
            velocityS *= model.getDensity();
            velocityS *= model.getVelocityS();

            temp = 2 * gradMu;
            temp *= model.getDensity();
            temp *= model.getVelocityS();

            velocityS += temp;
        }
    } else {
        this->initParameterisation(velocityS, ctx, dist, 0.0);
    }

    if (workflow.getInvertForDensity()) {
        if (model.getParameterisation() == 0) {
            density = gradRho;            
        } else if (model.getParameterisation() == 3) {
            density = scai::lama::pow(model.getVelocityP(), 2);
            temp = scai::lama::pow(model.getVelocityS(), 2);
            temp *= 2;
            density -= temp;
            density *= gradLambda;

            temp = scai::lama::pow(model.getVelocityS(), 2);
            temp *= gradMu;
            density += temp;

            density += gradRho;
        }
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }
          
    if (workflow.getInvertForPorosity() && (model.getParameterisation() == 1 || model.getParameterisation() == 2)) {            
        scai::lama::DenseVector<ValueType> rho_satDePorosity;
        scai::lama::DenseVector<ValueType> mu_satDePorosity;
        scai::lama::DenseVector<ValueType> K_satDePorosity; 
        // Based on Gassmann equation   
        // derivative of mu_sat with respect to porosity
        mu_satDePorosity = this->getMu_satDePorosity(model); 
        // sum with chain rule
        porosity = mu_satDePorosity * gradMu; 
        
        if (model.getParameterisation() == 2) { 
            // derivative of density with respect to porosity
            rho_satDePorosity = this->getDensityDePorosity(model);  
            // derivative of K_sat with respect to porosity
            K_satDePorosity = this->getK_satDePorosity(model); 
        
            // porosity derived from seismic modulus  
            rho_satDePorosity *= gradRho;  
            porosity += rho_satDePorosity;        
            K_satDePorosity *= gradK_sat;
            porosity += K_satDePorosity;  
        }
    } else {
        this->initParameterisation(porosity, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForSaturation() && model.getParameterisation() == 2) {        
        scai::lama::DenseVector<ValueType> rho_satDeSaturation;
        scai::lama::DenseVector<ValueType> K_satDeSaturation; 
        
        // Based on Gassmann equation     
        // derivative of density with respect to saturation
        rho_satDeSaturation = this->getDensityDeSaturation(model);   
        // derivative of K_sat with respect to saturation
        K_satDeSaturation = this->getK_satDeSaturation(model); 
        
        // water saturation derived from seismic modulus              
        // sum with chain rule
        saturation = rho_satDeSaturation * gradRho;       
        K_satDeSaturation *= gradK_sat;
        saturation += K_satDeSaturation; 
    } else {
        this->initParameterisation(saturation, ctx, dist, 0.0);
    }    
}

template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityStemp;
    scai::lama::DenseVector<ValueType> velocityPtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempX;
    scai::lama::DenseVector<ValueType> tempY;
    ValueType modelMax;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialFDorder = configEM.get<IndexType>("spatialFDorder");
        
    /* Get references to required derivatives matrixes */
    auto const &DxfEM = derivativesEM.getDxf();
    auto const &DyfEM = derivativesEM.getDyf();  
          
    densitytemp = model.getDensity(); 
    velocityStemp = model.getVelocityS();  
    velocityPtemp = model.getVelocityP();                  
                      
    scai::hmemo::ContextPtr ctx = densitytemp.getContextPtr();
    scai::dmemo::DistributionPtr dist = densitytemp.getDistributionPtr();  
          
    if (workflow.getInvertForDensity()) {
        densitytemp = modelTaper2DJoint.applyGradientTransformToEM(densitytemp);          
        
        tempX = DxfEM * densitytemp;    
        tempY = DyfEM * densitytemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder);   
        
        modelMax = tempX.maxNorm();
        if (modelMax != 0)
            tempX *= 1 / modelMax;  
        modelMax = tempY.maxNorm();
        if (modelMax != 0)
            tempY *= 1 / modelMax; 
        
        tempX = modelTaper2DJoint.applyGradientTransformToSeismic(tempX);  
        tempY = modelTaper2DJoint.applyGradientTransformToSeismic(tempY);                           
    } else {
        this->initParameterisation(tempX, ctx, dist, 0.0);
        this->initParameterisation(tempY, ctx, dist, 0.0);
    }   
    
    modelDerivativeXtemp = tempX;
    modelDerivativeYtemp = tempY;   
               
    if (workflow.getInvertForVs()) {
        velocityStemp = modelTaper2DJoint.applyGradientTransformToEM(velocityStemp); 
        
        tempX = DxfEM * velocityStemp;    
        tempY = DyfEM * velocityStemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder);  
        
        modelMax = tempX.maxNorm();
        if (modelMax != 0)
            tempX *= 1 / modelMax;  
        modelMax = tempY.maxNorm();
        if (modelMax != 0)
            tempY *= 1 / modelMax; 
        
        tempX = modelTaper2DJoint.applyGradientTransformToSeismic(tempX);    
        tempY = modelTaper2DJoint.applyGradientTransformToSeismic(tempY);                       
    } else {
        this->initParameterisation(tempX, ctx, dist, 0.0);
        this->initParameterisation(tempY, ctx, dist, 0.0);
    }    
    
    modelDerivativeXtemp += tempX;
    modelDerivativeYtemp += tempY;    
    
    if (workflow.getInvertForVp()) {
        velocityPtemp = modelTaper2DJoint.applyGradientTransformToEM(velocityPtemp);  
        
        tempX = DxfEM * velocityPtemp;    
        tempY = DyfEM * velocityPtemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder);  
        
        modelMax = tempX.maxNorm();
        if (modelMax != 0)
            tempX *= 1 / modelMax;  
        modelMax = tempY.maxNorm();
        if (modelMax != 0)
            tempY *= 1 / modelMax; 
        
        tempX = modelTaper2DJoint.applyGradientTransformToSeismic(tempX);  
        tempY = modelTaper2DJoint.applyGradientTransformToSeismic(tempY);   
    } else {
        this->initParameterisation(tempX, ctx, dist, 0.0);
        this->initParameterisation(tempY, ctx, dist, 0.0);
    }    
    
    modelDerivativeXtemp += tempX;
    modelDerivativeYtemp += tempY;    
    
    dataMisfit.setModelDerivativeX(modelDerivativeXtemp);
    dataMisfit.setModelDerivativeY(modelDerivativeYtemp);
}


template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityStemp;
    scai::lama::DenseVector<ValueType> velocityPtemp;
    scai::lama::DenseVector<ValueType> modelEMDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelEMDerivativeYtemp;
    scai::lama::DenseVector<ValueType> temp;
    scai::lama::DenseVector<ValueType> tempX;
    scai::lama::DenseVector<ValueType> tempY;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialFDorder = configEM.get<IndexType>("spatialFDorder");
    
    /* Get references to required derivatives matrixes */
    auto const &DxfEM = derivativesEM.getDxf();
    auto const &DyfEM = derivativesEM.getDyf();  
    
    densitytemp = model.getDensity();    
    velocityStemp = model.getVelocityS();  
    velocityPtemp = model.getVelocityP();  
    
    modelEMDerivativeXtemp = dataMisfitEM.getModelDerivativeX();  
    modelEMDerivativeYtemp = dataMisfitEM.getModelDerivativeY();          
                      
    scai::hmemo::ContextPtr ctx = densitytemp.getContextPtr();
    scai::dmemo::DistributionPtr dist = densitytemp.getDistributionPtr();  
    scai::dmemo::DistributionPtr distEM = modelEMDerivativeXtemp.getDistributionPtr();  
                
    if (workflow.getInvertForDensity()) {
        densitytemp = modelTaper2DJoint.applyGradientTransformToEM(densitytemp); 
        
        tempX = DxfEM * densitytemp;    
        tempY = DyfEM * densitytemp;          
        
        // cross-gradient of density and modelEMDerivative   
        densitytemp = modelEMDerivativeXtemp * tempY;   
        temp = modelEMDerivativeYtemp * tempX;   
        densitytemp -= temp;         
        
        KITGPI::Common::applyMedianFilterTo2DVector(densitytemp, NX, NY, spatialFDorder);
        
        density = modelTaper2DJoint.applyGradientTransformToSeismic(densitytemp);
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }    
                 
    if (workflow.getInvertForVs()) {
        velocityStemp = modelTaper2DJoint.applyGradientTransformToEM(velocityStemp); 
        
        tempX = DxfEM * velocityStemp;    
        tempY = DyfEM * velocityStemp;          
        
        // cross-gradient of vs and modelEMDerivative   
        velocityStemp = modelEMDerivativeXtemp * tempY;   
        temp = modelEMDerivativeYtemp * tempX;   
        velocityStemp -= temp;   
        
        KITGPI::Common::applyMedianFilterTo2DVector(velocityStemp, NX, NY, spatialFDorder);
        
        velocityS = modelTaper2DJoint.applyGradientTransformToSeismic(velocityStemp);      
    } else {
        this->initParameterisation(velocityS, ctx, dist, 0.0);
    }       
    
    if (workflow.getInvertForVp()) {
        velocityPtemp = modelTaper2DJoint.applyGradientTransformToEM(velocityPtemp);    
        
        tempX = DxfEM * velocityPtemp;    
        tempY = DyfEM * velocityPtemp;          
        
        // cross-gradient of vs and modelEMDerivative   
        velocityPtemp = modelEMDerivativeXtemp * tempY;   
        temp = modelEMDerivativeYtemp * tempX;   
        velocityPtemp -= temp;   
        
        KITGPI::Common::applyMedianFilterTo2DVector(velocityPtemp, NX, NY, spatialFDorder);
        
        velocityP = modelTaper2DJoint.applyGradientTransformToSeismic(velocityPtemp);          
    } else {
        this->initParameterisation(velocityP, ctx, dist, 0.0);
    }      
}

template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityStemp;
    scai::lama::DenseVector<ValueType> velocityPtemp;
    scai::lama::DenseVector<ValueType> modelEMDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelEMDerivativeYtemp;
    scai::lama::DenseVector<ValueType> temp;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialFDorder = configEM.get<IndexType>("spatialFDorder");
        
    densitytemp = this->getDensity();
    velocityStemp = this->getVelocityS();   
    velocityPtemp = this->getVelocityP();    
        
    modelEMDerivativeXtemp = dataMisfitEM.getModelDerivativeX();  
    modelEMDerivativeYtemp = dataMisfitEM.getModelDerivativeY();          
                      
    scai::hmemo::ContextPtr ctx = densitytemp.getContextPtr();
    scai::dmemo::DistributionPtr dist = densitytemp.getDistributionPtr();  
    scai::dmemo::DistributionPtr distEM = modelEMDerivativeXtemp.getDistributionPtr();  
                  
    temp = modelEMDerivativeYtemp - modelEMDerivativeXtemp; 
    
    if (workflow.getInvertForDensity()) {
        densitytemp = modelTaper2DJoint.applyGradientTransformToEM(densitytemp);  
        
        // derivative of cross-gradient with respect to density   
        densitytemp *= temp;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(densitytemp, NX, NY, spatialFDorder);
                
        density = modelTaper2DJoint.applyGradientTransformToSeismic(densitytemp); 
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }    
                 
    if (workflow.getInvertForVs()) {
        velocityStemp = modelTaper2DJoint.applyGradientTransformToEM(velocityStemp);          
        
        // derivative of cross-gradient with respect to vs  
        velocityStemp *= temp;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(velocityStemp, NX, NY, spatialFDorder);
        
        velocityS = modelTaper2DJoint.applyGradientTransformToSeismic(velocityStemp);       
    } else {
        this->initParameterisation(velocityS, ctx, dist, 0.0);
    }    
                 
    if (workflow.getInvertForVp()) { 
        velocityPtemp = modelTaper2DJoint.applyGradientTransformToEM(velocityPtemp);         
        
        // derivative of cross-gradient with respect to vs  
        velocityPtemp *= temp;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(velocityPtemp, NX, NY, spatialFDorder);
          
        velocityP = modelTaper2DJoint.applyGradientTransformToSeismic(velocityPtemp);          
    } else {
        this->initParameterisation(velocityP, ctx, dist, 0.0);
    }       
}

template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPriori, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{    
    scai::lama::DenseVector<ValueType> velocityPPrioritemp;
    scai::lama::DenseVector<ValueType> velocitySPrioritemp;
    scai::lama::DenseVector<ValueType> densityPrioritemp;
    
    velocityP = model.getVelocityP();  
    velocityS = model.getVelocityS();  
    density = model.getDensity();  
    velocityPPrioritemp = modelPriori.getVelocityP(); 
    velocitySPrioritemp = modelPriori.getVelocityS();  
    densityPrioritemp = modelPriori.getDensity();  
            
    scai::hmemo::ContextPtr ctx = velocityP.getContextPtr();
    scai::dmemo::DistributionPtr dist = velocityP.getDistributionPtr();
    
    if (workflow.getInvertForVp()) {
        velocityPPrioritemp = velocityP - velocityPPrioritemp;
        velocityP = this->calcStabilizingFunctionalGradientPerModel(velocityPPrioritemp, config, dataMisfit);
    } else {
        this->initParameterisation(velocityP, ctx, dist, 0.0);
    }    
    if (workflow.getInvertForVs()) {
        velocitySPrioritemp = velocityS - velocitySPrioritemp;
        velocityS = this->calcStabilizingFunctionalGradientPerModel(velocitySPrioritemp, config, dataMisfit);
    } else {
        this->initParameterisation(velocityS, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForDensity()) { 
        densityPrioritemp = density - densityPrioritemp;
        density = this->calcStabilizingFunctionalGradientPerModel(densityPrioritemp, config, dataMisfit);
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }  
}

/*! \brief Apply a median filter to filter the extrame value of the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::Viscoelastic<ValueType>::applyMedianFilter(KITGPI::Configuration::Configuration config)
{
    scai::lama::DenseVector<ValueType> density_temp;
    scai::lama::DenseVector<ValueType> velocityS_temp;
    scai::lama::DenseVector<ValueType> velocityP_temp;
    scai::lama::DenseVector<ValueType> porosity_temp;
    scai::lama::DenseVector<ValueType> saturation_temp;
    
    density_temp = this->getDensity();
    velocityS_temp = this->getVelocityS();
    velocityP_temp = this->getVelocityP();
    porosity_temp = this->getPorosity();
    saturation_temp = this->getSaturation();
    
    scai::IndexType NZ = config.get<IndexType>("NZ");
    if (NZ == 1) {
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity_temp.size() / NY;
        scai::IndexType spatialFDorder = config.get<IndexType>("spatialFDorder");
                
        KITGPI::Common::applyMedianFilterTo2DVector(density_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(velocityS_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(velocityP_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(porosity_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(saturation_temp, NX, NY, spatialFDorder);
        
        this->setDensity(density_temp);
        this->setVelocityS(velocityS_temp);
        this->setVelocityP(velocityP_temp);
        this->setPorosity(porosity_temp);    
        this->setSaturation(saturation_temp);
    }
}

template class KITGPI::Gradient::Viscoelastic<float>;
template class KITGPI::Gradient::Viscoelastic<double>;
