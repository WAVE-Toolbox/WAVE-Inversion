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
void KITGPI::Gradient::Acoustic<ValueType>::sumGradientPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd)
{
    auto distBig = density.getDistributionPtr();
    auto dist = gradientPerShot.getDensity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> recoverMatrix = model.getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd));
    recoverMatrix.assignTranspose(recoverMatrix);
        
    scai::lama::DenseVector<ValueType> temp;
    scai::lama::DenseVector<ValueType> weightingVector;
    scai::IndexType NY = modelCoordinates.getNY();
    
    weightingVector = gradientPerShot.calcWeightingVector(gradientPerShot.getVelocityP(), NY, shotInd);    
    temp = weightingVector * gradientPerShot.getVelocityP();  
    temp = recoverMatrix * temp;
    velocityP += temp; 

    temp = weightingVector * gradientPerShot.getDensity();  
    temp = recoverMatrix * temp;
    density += temp; 
    
    temp = weightingVector * gradientPerShot.getPorosity();  
    temp = recoverMatrix * temp;
    porosity += temp; 
    
    temp = weightingVector * gradientPerShot.getSaturation();  
    temp = recoverMatrix * temp;
    saturation += temp; 
    
    temp = weightingVector * gradientPerShot.getReflectivity();  
    temp = recoverMatrix * temp;
    reflectivity += temp; 
}

/*! \brief calculate Gaussian kernel
\param model model
\param modelCoordinates coordinate class object of the model
\param FCmain main frequency
*/
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::calcGaussianKernel(scai::dmemo::CommunicatorPtr commAll, KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Configuration::Configuration config)
{
    scai::IndexType smoothGradient = config.getAndCatch("smoothGradient", 0);
    if (smoothGradient != 0) {
        double start_t = common::Walltime::get();
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity.size() / NY; // NX is different in stream configuration
        ValueType DH = config.get<ValueType>("DH");
        scai::lama::DenseVector<ValueType> velocity; 
        velocity = model.getVelocityP();
        ValueType velocityMean = velocity.sum() / velocity.size(); 
        
        KITGPI::Common::calcGaussianKernelFor2DVector(porosity, GaussianKernel, PX, PY, NX, NY, DH, velocityMean, config.get<ValueType>("CenterFrequencyCPML"), smoothGradient);
        
        double end_t = common::Walltime::get();
        HOST_PRINT(commAll, "\nCalculate Gaussian kernel with matrix size = " << (NX+PX-1)*(NY+PY-1) << " x " << (NX+PX-1)*(NY+PY-1) << " and kernel size = " << PX << " x " << PY << " (" << PX*DH << " m x " << PY*DH << " m) in " << end_t - start_t << " sec.\n");
    }
}

/*! \brief Smooth gradient by Gaussian window
\param modelCoordinates coordinate class object of the model
*/
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::smooth(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config)
{
    if (config.get<IndexType>("NZ") == 1) {
        scai::IndexType smoothGradient = config.getAndCatch("smoothGradient", 0);
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity.size() / NY; // NX is different in stream configuration
        if (smoothGradient != 0) {
            HOST_PRINT(commAll, "\nApply Gaussian filter to gradient\n");
            scai::lama::DenseVector<ValueType> vector2Dpadded;
            if (workflowInner.getInvertForVp()) {
                KITGPI::Common::pad2DVector(velocityP, vector2Dpadded, NX, NY, PX, PY); 
                velocityP = GaussianKernel * vector2Dpadded;
            }
            if (workflowInner.getInvertForDensity()) {
                KITGPI::Common::pad2DVector(density, vector2Dpadded, NX, NY, PX, PY); 
                density = GaussianKernel * vector2Dpadded;
            }
            if (workflowInner.getInvertForPorosity()) {
                KITGPI::Common::pad2DVector(porosity, vector2Dpadded, NX, NY, PX, PY); 
                porosity = GaussianKernel * vector2Dpadded;
            }
            if (workflowInner.getInvertForSaturation()) {
                KITGPI::Common::pad2DVector(saturation, vector2Dpadded, NX, NY, PX, PY); 
                saturation = GaussianKernel * vector2Dpadded;
            }
            if (workflowInner.getInvertForReflectivity()) {
                KITGPI::Common::pad2DVector(reflectivity, vector2Dpadded, NX, NY, PX, PY); 
                reflectivity = GaussianKernel * vector2Dpadded;
            }
        }
    }
}

/*! \brief Smooth gradient by Gaussian window
\param modelCoordinates coordinate class object of the model
*/
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::applyMedianFilter(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config)
{
    if (config.get<IndexType>("NZ") == 1) {
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity.size() / NY; // NX is different in stream configuration
    
        HOST_PRINT(commAll, "\nApply median filter to gradient\n");
        scai::IndexType spatialLength = config.get<IndexType>("spatialFDorder");
        
        if (workflowInner.getInvertForVp())
            KITGPI::Common::applyMedianFilterTo2DVector(velocityP, NX, NY, spatialLength);
        if (workflowInner.getInvertForDensity())
            KITGPI::Common::applyMedianFilterTo2DVector(density, NX, NY, spatialLength);
        if (workflowInner.getInvertForPorosity())
            KITGPI::Common::applyMedianFilterTo2DVector(porosity, NX, NY, spatialLength);
        if (workflowInner.getInvertForSaturation())
            KITGPI::Common::applyMedianFilterTo2DVector(saturation, NX, NY, spatialLength);
        if (workflowInner.getInvertForReflectivity())
            KITGPI::Common::applyMedianFilterTo2DVector(reflectivity, NX, NY, spatialLength);
    }
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
                    mu = scai::lama::pow(model.getVelocityP(), 2);
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
                maxValue = config.getAndCatch("upperReflectivityTh", 1.0) - config.getAndCatch("lowerReflectivityTh", -1.0);
            }      
            reflectivity *= 1 / reflectivity.maxNorm() * maxValue;      
        }
    }
}

/*! \brief Function for applying EnergyPreconditioning to the gradient
 \param epsilonHessian a small value to stabilize the calculation of inverse Hessian.
 \param saveApproxHessian save the approximated Hessian.
 \param filename filename.
 \param fileFormat fileFormat.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::applyEnergyPreconditioning(ValueType epsilonHessian, scai::IndexType saveApproxHessian, std::string filename, scai::IndexType fileFormat)
{
    // see Nuber et al., 2015., Enhancement of near-surface elastic full waveform inversion results in regions of low sensitivities.
    scai::lama::DenseVector<ValueType> approxHessian;  
       
    if (workflowInner.getInvertForVs() && velocityP.maxNorm() != 0) {
        approxHessian = velocityP;
        approxHessian *= velocityP;
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, filename + ".vs", fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        velocityP *= approxHessian; 
    }
    
    if (workflowInner.getInvertForDensity() && density.maxNorm() != 0) {
        approxHessian = density;
        approxHessian *= density;
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, filename + ".density", fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        density *= approxHessian; 
    }
    
    if (workflowInner.getInvertForPorosity() && porosity.maxNorm() != 0) {
        approxHessian = porosity;
        approxHessian *= porosity;
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, filename + ".porosity", fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        porosity *= approxHessian; 
    }
    
    if (workflowInner.getInvertForSaturation() && saturation.maxNorm() != 0) {
        approxHessian = saturation;
        approxHessian *= saturation;
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, filename + ".saturation", fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        saturation *= approxHessian; 
    }
    
    if (workflowInner.getInvertForReflectivity() && reflectivity.maxNorm() != 0) {
        approxHessian = reflectivity;
        approxHessian *= reflectivity;
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, filename + ".reflectivity", fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        reflectivity *= approxHessian; 
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

/*! \brief Function for calculating the model derivatives 
 *
 \param dataMisfit dataMisfit to store the model derivatives.
 \param model model parameters.
 \param derivativesEM EM derivatives
 \param configEM EM config handle
 \param modelTaper2DJoint Taper to provide transform matrix
 \param workflow workflow
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityPtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempX;
    scai::lama::DenseVector<ValueType> tempY;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialLength = configEM.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = configEM.get<IndexType>("exchangeStrategy");
    ValueType DT = configEM.get<ValueType>("DT");
    
    /* Get references to required derivatives matrices */
    scai::lama::CSRSparseMatrix<ValueType> DxfEM;
    scai::lama::CSRSparseMatrix<ValueType> DyfEM;
    DxfEM = derivativesEM.getDxf();
    DyfEM = derivativesEM.getDyf();
    DxfEM.scale(1.0 / DT);
    DyfEM.scale(1.0 / DT);            
                      
    scai::hmemo::ContextPtr ctx = model.getDensity().getContextPtr();
    scai::dmemo::DistributionPtr dist = model.getDensity().getDistributionPtr();   
               
    densitytemp = modelTaper2DJoint.applyGradientTransformToEM(model.getDensity());
    velocityPtemp = modelTaper2DJoint.applyGradientTransformToEM(model.getVelocityP()); 
    
    // store the mean value of model parameters for weighting the gradient in summing
    if (workflow.workflowStage + workflow.iteration == 0) {
        density0mean = densitytemp.sum() / densitytemp.size();
        velocityP0mean = velocityPtemp.sum() / velocityPtemp.size();
    }
    
    velocityPtemp *= 1 / velocityP0mean;
    
    tempX = DxfEM * velocityPtemp;    
    tempY = DyfEM * velocityPtemp;
    
    KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialLength);
    KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialLength);  
        
    modelDerivativeXtemp = modelTaper2DJoint.applyGradientTransformToSeismic(tempX);    
    modelDerivativeYtemp = modelTaper2DJoint.applyGradientTransformToSeismic(tempY);  
        
    if (exchangeStrategy == 2) {   
        densitytemp *= 1 / density0mean;
        
        tempX = DxfEM * densitytemp;    
        tempY = DyfEM * densitytemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialLength);   
        
        modelDerivativeXtemp += modelTaper2DJoint.applyGradientTransformToSeismic(tempX);  
        modelDerivativeYtemp += modelTaper2DJoint.applyGradientTransformToSeismic(tempY); 
    }
    
    dataMisfit.setModelDerivativeX(modelDerivativeXtemp);
    dataMisfit.setModelDerivativeY(modelDerivativeYtemp);
}

/*! \brief Function for calculating the cross gradient 
 *
 \param dataMisfitEM EM dataMisfit to store the modelEM derivatives.
 \param model model parameters.
 \param derivativesEM EM derivatives
 \param configEM EM config handle
 \param modelTaper2DJoint Taper to provide transform matrix
 \param workflow workflow
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityPtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempXY;
    scai::lama::DenseVector<ValueType> tempYX;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialLength = configEM.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = configEM.get<IndexType>("exchangeStrategy");
    ValueType DT = configEM.get<ValueType>("DT");
    
    /* Get references to required derivatives matrices */
    scai::lama::CSRSparseMatrix<ValueType> DxfEM;
    scai::lama::CSRSparseMatrix<ValueType> DyfEM;
    DxfEM = derivativesEM.getDxf();
    DyfEM = derivativesEM.getDyf();
    DxfEM.scale(1.0 / DT);
    DyfEM.scale(1.0 / DT);
                      
    scai::hmemo::ContextPtr ctx = model.getDensity().getContextPtr();
    scai::dmemo::DistributionPtr dist = model.getDensity().getDistributionPtr();  
        
    densitytemp = modelTaper2DJoint.applyGradientTransformToEM(model.getDensity());  
    velocityPtemp = modelTaper2DJoint.applyGradientTransformToEM(model.getVelocityP()); 
        
    // store the mean value of model parameters for weighting the gradient in summing
    if (workflow.workflowStage + workflow.iteration == 0) {
        density0mean = densitytemp.sum() / densitytemp.size();
        velocityP0mean = velocityPtemp.sum() / velocityPtemp.size();
    }
    
    velocityPtemp *= 1 / velocityP0mean;
    
    modelDerivativeXtemp = DxfEM * velocityPtemp;    
    modelDerivativeYtemp = DyfEM * velocityPtemp;
    
    KITGPI::Common::applyMedianFilterTo2DVector(modelDerivativeXtemp, NX, NY, spatialLength);
    KITGPI::Common::applyMedianFilterTo2DVector(modelDerivativeYtemp, NX, NY, spatialLength);       
                                                   
    if (workflow.getInvertForVs() && exchangeStrategy != 0) {  
        // cross gradient of vs and EM model   
        tempYX = modelDerivativeYtemp * dataMisfitEM.getModelDerivativeX();   
        tempXY = modelDerivativeXtemp * dataMisfitEM.getModelDerivativeY();
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);   
        
        velocityPtemp = tempYX - tempXY;   
        KITGPI::Common::applyMedianFilterTo2DVector(velocityPtemp, NX, NY, spatialLength);
        
        velocityP = modelTaper2DJoint.applyGradientTransformToSeismic(velocityPtemp);  
    } else {
        this->initParameterisation(velocityP, ctx, dist, 0.0);
    }       
    
    if (workflow.getInvertForDensity()) { 
        // cross gradient of density and vs   
        densitytemp *= 1 / density0mean;
        
        tempXY = DxfEM * densitytemp;    
        tempYX = DyfEM * densitytemp; 
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);          
        
        tempYX = tempYX * modelDerivativeXtemp;   
        tempXY = tempXY * modelDerivativeYtemp;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);   
        
        densitytemp = tempYX - tempXY;    
        KITGPI::Common::applyMedianFilterTo2DVector(densitytemp, NX, NY, spatialLength);     
               
        density = modelTaper2DJoint.applyGradientTransformToSeismic(densitytemp);
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }        
}

/*! \brief Function for calculating the derivatives of cross gradient 
 *
 \param dataMisfitEM EM dataMisfit to store the modelEM derivatives.
 \param model model parameters.
 \param derivativesEM EM derivatives
 \param configEM EM config handle
 \param modelTaper2DJoint Taper to provide transform matrix
 \param workflow workflow
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityPtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempXY;
    scai::lama::DenseVector<ValueType> tempYX;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialLength = configEM.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = configEM.get<IndexType>("exchangeStrategy");
    ValueType DT = configEM.get<ValueType>("DT");
        
    scai::hmemo::ContextPtr ctx = model.getDensity().getContextPtr();
    scai::dmemo::DistributionPtr dist = model.getDensity().getDistributionPtr();  
                  
    /* Get references to required derivatives matrices */
    scai::lama::CSRSparseMatrix<ValueType> DxfEM;
    scai::lama::CSRSparseMatrix<ValueType> DyfEM;
    DxfEM = derivativesEM.getDxf();
    DyfEM = derivativesEM.getDyf();
    DxfEM.scale(1.0 / DT);
    DyfEM.scale(1.0 / DT);
                     
    if (workflow.getInvertForVs() && exchangeStrategy != 0) {    
        velocityPtemp = modelTaper2DJoint.applyGradientTransformToEM(velocityP);          
        
        // derivative of cross gradient of vs and EM model with respect to vs  
        tempYX = dataMisfitEM.getModelDerivativeX() * velocityPtemp; 
        tempXY = dataMisfitEM.getModelDerivativeY() * velocityPtemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);   
        
        tempYX = DyfEM * tempYX;   
        tempXY = DxfEM * tempXY;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);
        
        velocityPtemp = tempXY - tempYX;   
        KITGPI::Common::applyMedianFilterTo2DVector(velocityPtemp, NX, NY, spatialLength);
                
        velocityP = modelTaper2DJoint.applyGradientTransformToSeismic(velocityPtemp);    
    } else {
        this->initParameterisation(velocityP, ctx, dist, 0.0);
    }   
        
    velocityPtemp = modelTaper2DJoint.applyGradientTransformToEM(model.getVelocityP()); 
    
    velocityPtemp *= 1 / velocityP0mean;
    
    modelDerivativeXtemp = DxfEM * velocityPtemp;    
    modelDerivativeYtemp = DyfEM * velocityPtemp;
    
    KITGPI::Common::applyMedianFilterTo2DVector(modelDerivativeXtemp, NX, NY, spatialLength);
    KITGPI::Common::applyMedianFilterTo2DVector(modelDerivativeYtemp, NX, NY, spatialLength);       
                                   
    if (workflow.getInvertForDensity()) {
        densitytemp = modelTaper2DJoint.applyGradientTransformToEM(density);  
        
        // derivative of cross gradient of density and vs with respect to density  
        tempYX = modelDerivativeXtemp * densitytemp; 
        tempXY = modelDerivativeYtemp * densitytemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);
        
        tempYX = DyfEM * tempYX;   
        tempXY = DxfEM * tempXY;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);
        
        densitytemp = tempXY - tempYX;  
        KITGPI::Common::applyMedianFilterTo2DVector(densitytemp, NX, NY, spatialLength);
                       
        density = modelTaper2DJoint.applyGradientTransformToSeismic(densitytemp); 
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }         
}

/*! \brief calculate the misfit of CrossGradient
 *
 */
template <typename ValueType>
ValueType KITGPI::Gradient::Acoustic<ValueType>::calcCrossGradientMisfit()
{
    ValueType misfitSum = 0;
    IndexType count = 0;
    if (velocityP.l2Norm() != 0) {
        misfitSum += velocityP.l2Norm() / velocityP.size();
        count++;
    }
    if (density.l2Norm() != 0) {
        misfitSum += density.l2Norm() / density.size();
        count++;
    }
    if (count != 0)
        misfitSum /= count;
    
    return (misfitSum);
}

/*! \brief Function for calculating the stabilizing functional gradient 
 *
 \param model model parameters.
 \param modelPriori Priori model parameters.
 \param dataMisfit dataMisfit.
 \param workflow workflow.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPriori, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{    
    scai::lama::DenseVector<ValueType> velocityPPrioritemp;
    scai::lama::DenseVector<ValueType> densityPrioritemp;
    
    velocityP = model.getVelocityP();  
    density = model.getDensity();  
    velocityPPrioritemp = modelPriori.getVelocityP();  
    densityPrioritemp = modelPriori.getDensity();  
            
    scai::hmemo::ContextPtr ctx = velocityP.getContextPtr();
    scai::dmemo::DistributionPtr dist = velocityP.getDistributionPtr();
    
    if (workflow.getInvertForVs()) {
        velocityPPrioritemp = velocityP - velocityPPrioritemp;
        velocityP = this->calcStabilizingFunctionalGradientPerModel(velocityPPrioritemp, config, dataMisfit);
    } else {
        this->initParameterisation(velocityP, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForDensity()) { 
        densityPrioritemp = density - densityPrioritemp;
        density = this->calcStabilizingFunctionalGradientPerModel(densityPrioritemp, config, dataMisfit);
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }    
}

template class KITGPI::Gradient::Acoustic<double>;
template class KITGPI::Gradient::Acoustic<float>;
