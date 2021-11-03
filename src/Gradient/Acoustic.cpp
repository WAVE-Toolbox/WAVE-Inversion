#include "Acoustic.hpp"
#include <scai/dmemo/SingleDistribution.hpp>

using namespace scai;
using namespace KITGPI;

/*! \brief Constructor that is using the Configuration class
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType>::Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "acoustic";
    init(ctx, dist, 0.0, 0.0);
}

/*! \brief Initialisation that is using the configuration class
 *
 *  Generates a homogeneous gradient, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist, 0.0, 0.0);
}

/*! \brief Initialisation that is generating a homogeneous gradient
 *
 *  Generates a homogeneous gradient, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param velocityP_const velocity gradients given as scalar value
 \param rho_const Density gradient given as scalar value
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType rho_const)
{
    this->initParameterisation(velocityP, ctx, dist, velocityP_const);
    this->initParameterisation(density, ctx, dist, rho_const);
    this->initParameterisation(porosity, ctx, dist, 0.0);
    this->initParameterisation(saturation, ctx, dist, 0.0);
    this->initParameterisation(reflectivity, ctx, dist, 0.0);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType>::Acoustic(const Acoustic &rhs)
{
    equationType = rhs.equationType;
    velocityP = rhs.velocityP;
    density = rhs.density;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    reflectivity = rhs.reflectivity;
}

/*! \brief Set all parameter to zero.
*/
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::resetGradient()
{
    this->resetParameter(velocityP);
    this->resetParameter(density);
    this->resetParameter(porosity);
    this->resetParameter(saturation);
    this->resetParameter(reflectivity);
}
/*! \brief Write gradient to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added and for density ".density.mtx" is added.
 \param fileFormat format of output file
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::write(std::string filename, IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const
{
    if (workflow.getInvertForVp()) {
        std::string filenameP = filename + ".vp";
        this->writeParameterisation(velocityP, filenameP, fileFormat);
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
    
    if (workflow.getInvertForReflectivity()) { 
        std::string filenameReflectivity = filename + ".reflectivity";
        this->writeParameterisation(reflectivity, filenameReflectivity, fileFormat);
    }
};

/*! \brief Get equationType (acoustic)
 */
template <typename ValueType>
std::string KITGPI::Gradient::Acoustic<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Acoustic<ValueType>::getVelocityS()
{
    COMMON_THROWEXCEPTION("The S-wave velocity is not defined in an acoustic simulation.")
    return (velocityS);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Acoustic<ValueType>::getTauP()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Acoustic<ValueType>::getTauS()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauS);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Gradient::Acoustic<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an elastic modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Gradient::Acoustic<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an elastic modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Overloading * Operation
 *
 \param rhs scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> KITGPI::Gradient::Acoustic<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Gradient::Acoustic<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief non-member function to multiply (scalar as left operand)
 *
 \param lhs scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> operator*(ValueType lhs, KITGPI::Gradient::Acoustic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> &KITGPI::Gradient::Acoustic<ValueType>::operator*=(ValueType const &rhs)
{
    if (workflowInner.getInvertForDensity()) 
        density *= rhs;
    if (workflowInner.getInvertForVp()) 
        velocityP *= rhs;
    if (workflowInner.getInvertForPorosity()) 
        porosity *= rhs;
    if (workflowInner.getInvertForSaturation()) 
        saturation *= rhs;
    if (workflowInner.getInvertForReflectivity()) 
        reflectivity *= rhs;

    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> KITGPI::Gradient::Acoustic<ValueType>::operator+(KITGPI::Gradient::Acoustic<ValueType> const &rhs)
{
    KITGPI::Gradient::Acoustic<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> &KITGPI::Gradient::Acoustic<ValueType>::operator+=(KITGPI::Gradient::Acoustic<ValueType> const &rhs)
{
    density += rhs.density;
    velocityP += rhs.velocityP;
    porosity += rhs.porosity;
    saturation += rhs.saturation;
    reflectivity += rhs.reflectivity;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Gradient which is subtractet.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> KITGPI::Gradient::Acoustic<ValueType>::operator-(KITGPI::Gradient::Acoustic<ValueType> const &rhs)
{
    KITGPI::Gradient::Acoustic<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> &KITGPI::Gradient::Acoustic<ValueType>::operator-=(KITGPI::Gradient::Acoustic<ValueType> const &rhs)
{
    density -= rhs.density;
    velocityP -= rhs.velocityP;
    porosity -= rhs.porosity;
    saturation -= rhs.saturation;
    reflectivity -= rhs.reflectivity;
    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Gradient which is copied.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> &KITGPI::Gradient::Acoustic<ValueType>::operator=(KITGPI::Gradient::Acoustic<ValueType> const &rhs)
{
    // why does rhs.density not work (density = protected)
    density = rhs.density;
    velocityP = rhs.velocityP;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    reflectivity = rhs.reflectivity;
    return *this;
}

/*! \brief function for overloading = Operation (called in base class)
 *
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::assign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    density = rhs.getDensity();
    velocityP = rhs.getVelocityP();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    reflectivity = rhs.getReflectivity();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    density -= rhs.getDensity();
    velocityP -= rhs.getVelocityP();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    reflectivity -= rhs.getReflectivity();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    density += rhs.getDensity();
    velocityP += rhs.getVelocityP();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    reflectivity += rhs.getReflectivity();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::timesAssign(ValueType const &rhs)
{
    if (workflowInner.getInvertForDensity()) 
        density *= rhs;
    if (workflowInner.getInvertForVp()) 
        velocityP *= rhs;
    if (workflowInner.getInvertForPorosity()) 
        porosity *= rhs;
    if (workflowInner.getInvertForSaturation()) 
        saturation *= rhs;
    if (workflowInner.getInvertForReflectivity()) 
        reflectivity *= rhs;
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::timesAssign(scai::lama::Vector<ValueType> const &rhs)
{
    density *= rhs;
    velocityP *= rhs;
    porosity *= rhs;
    saturation *= rhs;
    reflectivity *= rhs;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract gradient.
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs)
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
        if (workflowInner.getInvertForDensity()) {
            temp = lhs.getDensity() - rhs.getDensity();
            lhs.setDensity(temp);
        }
    } else if (lhs.getParameterisation() == 0) {
        scai::lama::DenseVector<ValueType> rho;
        rho = lhs.getDensity();
        if (workflowInner.getInvertForDensity()) {
            temp = lhs.getDensity() - rhs.getDensity();
            lhs.setDensity(temp);
        }        
        if (workflowInner.getInvertForVp()) {
            scai::lama::DenseVector<ValueType> lambda;
            lambda = scai::lama::pow(lhs.getVelocityP(), 2);
            lambda *= rho;        
            lambda -= rhs.getVelocityP();  
            lambda /= rho;
            lambda = scai::lama::sqrt(lambda);
            Common::replaceInvalid<ValueType>(lambda, 0.0);
            lhs.setVelocityP(lambda);
        }
    }
    if (workflowInner.getInvertForReflectivity()) {
        temp = lhs.getReflectivity() - rhs.getReflectivity();
        lhs.setReflectivity(temp); 
    }
}

/*! \brief Function for summing the gradients of all shot domains
 *
 \param commInterShot inter shot communication pointer
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot)
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
    
    density.redistribute(singleDist);
    commInterShot->sumArray(density.getLocalValues());
    density.redistribute(dist);
    
    porosity.redistribute(singleDist);
    commInterShot->sumArray(porosity.getLocalValues());
    porosity.redistribute(dist);
    
    saturation.redistribute(singleDist);
    commInterShot->sumArray(saturation.getLocalValues());
    saturation.redistribute(dist);
    
    reflectivity.redistribute(singleDist);
    commInterShot->sumArray(reflectivity.getLocalValues());
    reflectivity.redistribute(dist);
}

/*! \brief If stream configuration is used, get a pershot model from the big model
 \param model model
 \param gradientPerShot gradient per shot
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::sumGradientPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth)
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
    
    temp = recoverMatrix * gradientPerShot.getReflectivity(); //transform pershot into big model
    temp *= recoverVector;
    reflectivity *= eraseVector;
    reflectivity += temp; //take over the values
}

/*! \brief function for scaling the gradients with the model parameter
 *
 \param model Abstract model.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{
    ValueType maxValue = 1;      
    
    IndexType scaleGradient = config.get<IndexType>("scaleGradient");
    if (scaleGradient != 0) {
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
                    lambdaMax = pow(config.get<ValueType>("upperVPTh"), 2);
                    lambdaMax *= config.get<ValueType>("upperDensityTh");
                    lambdaMin = pow(config.get<ValueType>("lowerVPTh"), 2);
                    lambdaMin *= config.get<ValueType>("lowerDensityTh");
                    maxValue = lambdaMax - lambdaMin;
                }
            }
            velocityP *= 1 / velocityP.maxNorm() * maxValue;
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
            
        if (workflow.getInvertForReflectivity() && reflectivity.maxNorm() != 0) { 
            if (scaleGradient == 1) {
                maxValue = model.getReflectivity().maxNorm();
            } else if (scaleGradient == 2) {
                maxValue = config.getAndCatch("upperReflectivityEMrTh", 1.0) - config.getAndCatch("lowerReflectivityEMrTh", -1.0);
            }      
            reflectivity *= 1 / reflectivity.maxNorm() * maxValue;      
        }
    }
}

/*! \brief Function for normalizing the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::normalize()
{
    if (this->getNormalizeGradient()) {
        ValueType gradientMax = velocityP.maxNorm();
        if (gradientMax != 0)
            velocityP *= 1 / gradientMax;
        gradientMax = density.maxNorm();
        if (gradientMax != 0)
            density *= 1 / gradientMax;
        if (gradientMax != 0)
            porosity *= 1 / gradientMax;
        gradientMax = saturation.maxNorm();
        if (gradientMax != 0)
            saturation *= 1 / gradientMax;
        gradientMax = reflectivity.maxNorm();
        if (gradientMax != 0)
            reflectivity *= 1 / gradientMax;
    }
}

/*! \brief function for calculating the acoustic gradients from the cross-correlation and the model parameter
 *
 \param model Abstract model.
 \param correlatedWavefields Abstract xCorr.
 \param DT Temporal discretization
 \param workflow
 *
 \f{eqnarray*}
   \nabla_K E &=& -\mathrm{d}t \frac{1}{v_{\mathrm{p}}^4 \rho^2} \cdot X_{\lambda} \\
   \nabla_{v_{\mathrm{p}}} E &=& 2 v_{\mathrm{p}} \rho \cdot \nabla_K E \\
   \nabla_{\rho} E &=& v_{\mathrm{p}}^2 \cdot \nabla_K E - \mathrm{d}t \cdot X_{\rho}
 \f}
 *
 \sa{KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::update} for the cross-correlations \f$ (X_{\lambda},X_{\rho}) \f$ in 2D
 \sa{KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::update} for the cross-correlations \f$ (X_{\lambda},X_{\rho}) \f$ in 3D
 *
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //dt should be in cross correlation!
    //gradBulk = -dt * Padj * dPfw/dt / lambda^2
    scai::lama::DenseVector<ValueType> gradBulk;
    scai::lama::DenseVector<ValueType> gradRho;
    
    gradBulk = scai::lama::pow(model.getVelocityP(), 2);
    gradBulk *= model.getDensity();
    gradBulk = scai::lama::pow(gradBulk, -2);
    Common::replaceInvalid<ValueType>(gradBulk, 0.0);

    gradBulk *= correlatedWavefields.getXcorrLambda();

    gradBulk *= -DT;
    
    scai::hmemo::ContextPtr ctx = gradBulk.getContextPtr();
    scai::dmemo::DistributionPtr dist = gradBulk.getDistributionPtr();

    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        gradRho = -DT * correlatedWavefields.getXcorrRho(); 
    }
    
    if (workflow.getInvertForVp()) {
        if (model.getParameterisation() == 0) {
            velocityP = gradBulk;
        } else if (model.getParameterisation() == 3) {
            //grad_vp = 2 * rho *vp * gradBulk
            velocityP = 2 * gradBulk;
            velocityP *= model.getDensity();
            velocityP *= model.getVelocityP();
        }
    } else {
        this->initParameterisation(velocityP, ctx, dist, 0.0);
    }

    if (workflow.getInvertForDensity()) {
        if (model.getParameterisation() == 0) {
            density = gradRho;            
        } else if (model.getParameterisation() == 3) {
            //grad_density' = vp^2 * gradBulk + (-) (sum(i) (dt Vadj,i * dVfw,i/dt))
            density = scai::lama::pow(model.getVelocityP(), 2);
            density *= gradBulk;
            density += gradRho;
        }
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForPorosity() && (model.getParameterisation() == 1 || model.getParameterisation() == 2)) {            
        scai::lama::DenseVector<ValueType> rho_satDePorosity;
        scai::lama::DenseVector<ValueType> K_satDePorosity; 
        // Based on Gassmann equation   
        // derivative of K_sat with respect to porosity
        K_satDePorosity = this->getK_satDePorosity(model); 
        // sum with chain rule
        porosity = K_satDePorosity * gradBulk; 
        
        if (model.getParameterisation() == 2) { 
            // derivative of density with respect to porosity
            rho_satDePorosity = this->getDensityDePorosity(model);  
        
            // porosity derived from seismic modulus  
            rho_satDePorosity *= gradRho;  
            porosity += rho_satDePorosity;   
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
        K_satDeSaturation *= gradBulk;
        saturation += K_satDeSaturation; 
    } else {
        this->initParameterisation(saturation, ctx, dist, 0.0);
    }   
    
    if (workflow.getInvertForReflectivity()) {
        reflectivity = correlatedWavefields.getXcorrLambda();
        reflectivity *= -1;
    } else {
        this->initParameterisation(reflectivity, ctx, dist, 0.0);
    }  
}

template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
}

template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    
}

template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
}

/*! \brief calculate the misfit of CrossGradient
 *
 */
template <typename ValueType>
ValueType KITGPI::Gradient::Acoustic<ValueType>::calcCrossGradientMisfit()
{
    ValueType misfitSum = velocityP.l2Norm() + density.l2Norm();
    
    return (misfitSum);
}

template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPriori, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{    

}

/*! \brief Apply a median filter to filter the extrame value of the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::applyMedianFilter(KITGPI::Configuration::Configuration config)
{
    scai::lama::DenseVector<ValueType> density_temp;
    scai::lama::DenseVector<ValueType> velocityP_temp;
    scai::lama::DenseVector<ValueType> porosity_temp;
    scai::lama::DenseVector<ValueType> saturation_temp;
    scai::lama::DenseVector<ValueType> reflectivity_temp;
    
    density_temp = this->getDensity();
    velocityP_temp = this->getVelocityP();
    porosity_temp = this->getPorosity();
    saturation_temp = this->getSaturation();
    reflectivity_temp = this->getReflectivity();

    scai::IndexType NZ = config.get<IndexType>("NZ");
    if (NZ == 1) {
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity_temp.size() / NY;
        scai::IndexType spatialFDorder = config.get<IndexType>("spatialFDorder");
            
        KITGPI::Common::applyMedianFilterTo2DVector(density_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(velocityP_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(porosity_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(saturation_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(reflectivity_temp, NX, NY, spatialFDorder);
        
        this->setDensity(density_temp);
        this->setVelocityP(velocityP_temp);
        this->setPorosity(porosity_temp);    
        this->setSaturation(saturation_temp);
        this->setReflectivity(reflectivity_temp);
    }
}

template class KITGPI::Gradient::Acoustic<double>;
template class KITGPI::Gradient::Acoustic<float>;
