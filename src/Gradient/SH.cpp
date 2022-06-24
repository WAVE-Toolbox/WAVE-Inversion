#include "SH.hpp"
#include <scai/dmemo/SingleDistribution.hpp>

using namespace scai;

/*! \brief Constructor that is using the Configuration class
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Gradient::SH<ValueType>::SH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "sh";
    init(ctx, dist, 0.0, 0.0);
}

/*! \brief Initialisation with zeros
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist, 0.0, 0.0);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param velocityS_const S-wave modulus given as Scalar
 \param rho_const Density given as Scalar
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityS_const, ValueType rho_const)
{
    this->initParameterisation(velocityS, ctx, dist, velocityS_const);
    this->initParameterisation(density, ctx, dist, rho_const);
    this->initParameterisation(porosity, ctx, dist, 0.0);
    this->initParameterisation(saturation, ctx, dist, 0.0);
    this->initParameterisation(reflectivity, ctx, dist, 0.0);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::SH<ValueType>::SH(const SH &rhs)
{
    equationType = rhs.equationType;
    velocityS = rhs.velocityS;
    density = rhs.density;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    reflectivity = rhs.reflectivity;
}

/*! \brief Set all parameter to zero.
*/
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::resetGradient()
{
    this->resetParameter(velocityS);
    this->resetParameter(density);
    this->resetParameter(porosity);
    this->resetParameter(saturation);
    this->resetParameter(reflectivity);
}

/*! \brief Write model to an external file
 *
 \param filename For the second ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 \param fileFormat format of output file
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::write(std::string filename, IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const
{
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
    
    if (workflow.getInvertForReflectivity()) { 
        std::string filenameReflectivity = filename + ".reflectivity";
        this->writeParameterisation(reflectivity, filenameReflectivity, fileFormat);
    }
};

/*! \brief Get equationType (sh)
 */
template <typename ValueType>
std::string KITGPI::Gradient::SH<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to velocityP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::SH<ValueType>::getVelocityP()
{
    COMMON_THROWEXCEPTION("There is no velocityP in an sh modelling")
    return (velocityP);
}
/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::SH<ValueType>::getTauP()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an sh modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::SH<ValueType>::getTauS()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an sh modelling")
    return (tauS);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Gradient::SH<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an sh modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Gradient::SH<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an sh modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::SH<ValueType> KITGPI::Gradient::SH<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Gradient::SH<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief Free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Gradient::SH<ValueType> operator*(ValueType lhs, KITGPI::Gradient::SH<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::SH<ValueType> &KITGPI::Gradient::SH<ValueType>::operator*=(ValueType const &rhs)
{
    if (workflowInner.getInvertForDensity()) 
        density *= rhs;
    if (workflowInner.getInvertForVs()) 
        velocityS *= rhs;
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
KITGPI::Gradient::SH<ValueType> KITGPI::Gradient::SH<ValueType>::operator+(KITGPI::Gradient::SH<ValueType> const &rhs)
{
    KITGPI::Gradient::SH<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::SH<ValueType> &KITGPI::Gradient::SH<ValueType>::operator+=(KITGPI::Gradient::SH<ValueType> const &rhs)
{
    density += rhs.density;
    velocityS += rhs.velocityS;
    porosity += rhs.porosity;
    saturation += rhs.saturation;
    reflectivity += rhs.reflectivity;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::SH<ValueType> KITGPI::Gradient::SH<ValueType>::operator-(KITGPI::Gradient::SH<ValueType> const &rhs)
{
    KITGPI::Gradient::SH<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::SH<ValueType> &KITGPI::Gradient::SH<ValueType>::operator-=(KITGPI::Gradient::SH<ValueType> const &rhs)
{
    density -= rhs.density;
    velocityS -= rhs.velocityS;
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
KITGPI::Gradient::SH<ValueType> &KITGPI::Gradient::SH<ValueType>::operator=(KITGPI::Gradient::SH<ValueType> const &rhs)
{
    density = rhs.density;
    velocityS = rhs.velocityS;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    reflectivity = rhs.reflectivity;

    return *this;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::assign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    density = rhs.getDensity();
    velocityS = rhs.getVelocityS();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    reflectivity = rhs.getReflectivity();
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    density -= rhs.getDensity();
    velocityS -= rhs.getVelocityS();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    reflectivity -= rhs.getReflectivity();
}

/*! \brief Function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    density += rhs.getDensity();
    velocityS += rhs.getVelocityS();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    reflectivity += rhs.getReflectivity();
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::timesAssign(ValueType const &rhs)
{
    if (workflowInner.getInvertForDensity()) 
        density *= rhs;
    if (workflowInner.getInvertForVs()) 
        velocityS *= rhs;
    if (workflowInner.getInvertForPorosity()) 
        porosity *= rhs;
    if (workflowInner.getInvertForSaturation()) 
        saturation *= rhs;
    if (workflowInner.getInvertForReflectivity()) 
        reflectivity *= rhs;
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::timesAssign(scai::lama::Vector<ValueType> const &rhs)
{
    density *= rhs;
    velocityS *= rhs;
    porosity *= rhs;
    saturation *= rhs;
    reflectivity *= rhs;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract model.
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs)
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
    }
    if (workflowInner.getInvertForReflectivity()) {
        temp = lhs.getReflectivity() - rhs.getReflectivity();
        lhs.setReflectivity(temp); 
    }
};

/*! \brief Function for summing the gradients of all shot domains
 *
 \param commInterShot inter shot communication pointer
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot)
{
    /*reduction between shot domains.
    each shot domain may have a different distribution of (gradient) vectors. This happens if geographer is used (different result for dist on each shot domain even for homogenous architecture) or on heterogenous architecture. In this case even the number of processes on each domain can vary. Therefore it is necessary that only one process per shot domain communicates all data.
    */
    
    //get information from distributed vector
    auto size = velocityS.size();
    auto dist = velocityS.getDistributionPtr();
    auto comm = dist->getCommunicatorPtr();
    
    // create single distribution, only master process owns the complete vector (no distribution).
    
    int shotMaster = 0;
    auto singleDist = std::make_shared<dmemo::SingleDistribution>( size, comm, shotMaster );
    
    //redistribute vector to master process
    // (this may cause memory issues for big models)
    velocityS.redistribute(singleDist);
    
    //reduce local array (size of local array is !=0 only for master process)
    commInterShot->sumArray(velocityS.getLocalValues());
    
    //redistribute vector to former partition
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
void KITGPI::Gradient::SH<ValueType>::sumGradientPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D cutCoordinate)
{
    auto distBig = density.getDistributionPtr();
    auto dist = gradientPerShot.getDensity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> recoverMatrix = model.getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);
    recoverMatrix.assignTranspose(recoverMatrix);
        
    scai::lama::DenseVector<ValueType> temp;
    scai::lama::DenseVector<ValueType> weightingVector;
    scai::IndexType NY = modelCoordinates.getNY();
    
    weightingVector = gradientPerShot.calcWeightingVector(gradientPerShot.getVelocityS(), NY);  
    
    temp = weightingVector * gradientPerShot.getVelocityS();  
    temp = recoverMatrix * temp;
    velocityS += temp; 

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

/*! \brief Smooth gradient by Gaussian window
\param modelCoordinates coordinate class object of the model
*/
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::smooth(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config)
{
    if (config.get<IndexType>("NZ") == 1) {
        scai::IndexType smoothGradient = config.getAndCatch("smoothGradient", 0);
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity.size() / NY; // NX is different in stream configuration
        if (smoothGradient != 0) {
            HOST_PRINT(commAll, "Apply Gaussian filter to gradient\n");
            scai::lama::DenseVector<ValueType> vector2Dpadded;
            if (workflowInner.getInvertForVs()) {
                KITGPI::Common::pad2DVector(velocityS, vector2Dpadded, NX, NY, PX, PY); 
                velocityS = GaussianKernel * vector2Dpadded;
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
void KITGPI::Gradient::SH<ValueType>::applyMedianFilter(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config)
{
    if (config.get<IndexType>("NZ") == 1) {
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity.size() / NY; // NX is different in stream configuration
    
        HOST_PRINT(commAll, "Apply median filter to gradient\n");
        scai::IndexType spatialLength = config.get<IndexType>("spatialFDorder");
        
        if (workflowInner.getInvertForVs())
            KITGPI::Common::applyMedianFilterTo2DVector(velocityS, NX, NY, spatialLength);
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

/*! \brief Function for scaling the gradients with the model parameter 
 *
 \param model Abstract model.
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{    
    ValueType maxValue = 1;      
    
    IndexType scaleGradient = config.get<IndexType>("scaleGradient");  
    if (scaleGradient != 0) {  
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
void KITGPI::Gradient::SH<ValueType>::applyEnergyPreconditioning(ValueType epsilonHessian, scai::IndexType saveApproxHessian, std::string filename, scai::IndexType fileFormat)
{
    // see Nuber et al., 2015., Enhancement of near-surface elastic full waveform inversion results in regions of low sensitivities.
    scai::lama::DenseVector<ValueType> approxHessian;  
       
    if (workflowInner.getInvertForVs() && velocityS.maxNorm() != 0) {
        approxHessian = velocityS;
        approxHessian *= velocityS;
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, filename + ".vs", fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        velocityS *= approxHessian; 
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
void KITGPI::Gradient::SH<ValueType>::normalize()
{
    if (this->getNormalizeGradient()) {
        ValueType gradientMax = velocityS.maxNorm();
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
        gradientMax = reflectivity.maxNorm();
        if (gradientMax != 0)
            reflectivity *= 1 / gradientMax;
    }
}

/*! \brief Function for calculating the sh gradients from the cross correlation and the model parameter 
 *
 \param model Abstract model.
 \param correlatedWavefields Abstract xCorr.
 \param DT Temporal discretization 
 \param workflow 
 *
 \f{eqnarray*}
   \nabla_{\mu} E &=& - \mathrm{d}t \frac{1}{(N \mu+2\mu)^2} \cdot X_{\mu} \\
   \nabla_{\mu} E &=& \mathrm{d}t \left[ \frac{N \mu^2 + 4 \mu \mu}{2 \mu^2 (N \mu+2\mu)^2} - \frac{1}{2 \mu^2} \right] \cdot X_{\mu,A} + \mathrm{d}t \frac{N \mu^2 + 4 \mu \mu}{2 \mu^2 (N \mu+2\mu)^2} \cdot X_{\mu,B} - \mathrm{d}t  \frac{1}{\mu^2} \cdot X_{\mu,C} \\ \\
   \nabla_{v_{\mathrm{p}}} E &=& 2 \rho v_{\mathrm{p}} \cdot \nabla_{\mu} E \\
   \nabla_{v_{\mathrm{s}}} E &=& 2 \rho v_{\mathrm{s}} \cdot \nabla_{\mu} E - 4 \rho v_{\mathrm{s}} \nabla_{\mu} E \\
   \nabla_{\rho} E &=& ( v_{\mathrm{p}}^2-2v_{\mathrm{s}}^2 ) \cdot \nabla_{\mu} E + v_{\mathrm{s}}^2 \cdot \nabla_{\mu} E - \mathrm{d}t \cdot X_{\rho}
 \f}
 *
 * with \f$ N \f$ as the number of dimensions, \f$ \mu = \rho v_{\mathrm{s}}^2 \f$ and \f$ \mu = \rho v_{\mathrm{p}}^2 - 2\mu\f$.
 *
 \sa{KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dsh<ValueType>::update} for the cross-correlations \f$ (X_{\mu},X_{\rho},X_{\mu,A},X_{\mu,B},X_{\mu,C}) \f$ in 2D 
 \sa{KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dsh<ValueType>::update} for the cross-correlations \f$ (X_{\mu},X_{\rho},X_{\mu,A},X_{\mu,B},X_{\mu,C}) \f$ in 3D
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //dt should be in cross correlation!
    //gradMu = -dt * Padj * dPfw/dt / mu^2
    scai::lama::DenseVector<ValueType> gradMu;
    scai::lama::DenseVector<ValueType> gradRho;
    
    gradMu = scai::lama::pow(model.getVelocityS(), 2);
    gradMu *= model.getDensity();
    gradMu = scai::lama::pow(gradMu, -2);
    Common::replaceInvalid<ValueType>(gradMu, 0.0);

    gradMu *= correlatedWavefields.getXcorrMuC();
    gradMu *= -DT;        
        
    scai::hmemo::ContextPtr ctx = gradMu.getContextPtr();
    scai::dmemo::DistributionPtr dist = gradMu.getDistributionPtr();
             
    if (workflow.getInvertForDensity() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        gradRho = -DT * correlatedWavefields.getXcorrRho(); 
    }
    
    if (workflow.getInvertForVs()) {
        if (model.getParameterisation() == 0) {
            velocityS = gradMu;
        } else if (model.getParameterisation() == 3) {
            //grad_vs = 2 * rho *vs * gradMu
            velocityS = 2 * gradMu;
            velocityS *= model.getDensity();
            velocityS *= model.getVelocityS();
        }
    } else {
        this->initParameterisation(velocityS, ctx, dist, 0.0);
    }

    if (workflow.getInvertForDensity()) {
        if (model.getParameterisation() == 0) {
            density = gradRho;            
        } else if (model.getParameterisation() == 3) {
            //grad_density' = vs^2 * gradMu + (-) (sum(i) (dt Vadj,i * dVfw,i/dt))
            density = scai::lama::pow(model.getVelocityS(), 2);
            density *= gradMu;
            density += gradRho;
        }
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForPorosity() && (model.getParameterisation() == 1 || model.getParameterisation() == 2)) {       
        scai::lama::DenseVector<ValueType> rho_satDePorosity;
        scai::lama::DenseVector<ValueType> mu_satDePorosity;
        // Based on Gassmann equation   
        // derivative of mu_sat with respect to porosity
        mu_satDePorosity = this->getMu_satDePorosity(model); 
        // sum with chain rule
        porosity = mu_satDePorosity * gradMu; 
        
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
        
        // Based on Gassmann equation     
        // derivative of density with respect to saturation
        rho_satDeSaturation = this->getDensityDeSaturation(model);  
        
        // water saturation derived from seismic modulus              
        // sum with chain rule
        saturation = rho_satDeSaturation * gradRho;   
    } else {
        this->initParameterisation(saturation, ctx, dist, 0.0);
    }    
    
    if (workflow.getInvertForReflectivity()) {
        reflectivity = correlatedWavefields.getXcorrMuC();
        reflectivity *= -1;
    } else {
        this->initParameterisation(reflectivity, ctx, dist, 0.0);
    }  
}

/*! \brief Function for calculating the model derivatives 
 *
 \param dataMisfit dataMisfit to store the model derivatives.
 \param model model parameters.
 \param derivatives derivatives
 \param config config handle
 \param modelTaper2DJoint Taper to provide transform matrix
 \param workflow workflow
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, KITGPI::Configuration::Configuration config, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityStemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempX;
    scai::lama::DenseVector<ValueType> tempY;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(config, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(config, "NY");
    scai::IndexType spatialLength = config.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = config.get<IndexType>("exchangeStrategy");
    
    /* Get references to required derivatives matrices */
    scai::lama::CSRSparseMatrix<ValueType> Dxf;
    scai::lama::CSRSparseMatrix<ValueType> Dyf;
    ValueType DT = config.get<ValueType>("DT");
    Dxf = derivatives.getDxf();
    Dyf = derivatives.getDyf();
    Dxf.scale(1.0 / DT);
    Dyf.scale(1.0 / DT);
                      
    densitytemp = model.getDensity();
    velocityStemp = model.getVelocityS(); 
    
    // store the mean value of model parameters for weighting the gradient in summing
    if (workflow.workflowStage + workflow.iteration == 0) {
        density0mean = densitytemp.sum() / densitytemp.size();
        velocityS0mean = velocityStemp.sum() / velocityStemp.size();
    }
    
    velocityStemp *= 1 / velocityS0mean;
    
    modelDerivativeXtemp = Dxf * velocityStemp;    
    modelDerivativeYtemp = Dyf * velocityStemp;
    
    KITGPI::Common::applyMedianFilterTo2DVector(modelDerivativeXtemp, NX, NY, spatialLength);
    KITGPI::Common::applyMedianFilterTo2DVector(modelDerivativeYtemp, NX, NY, spatialLength);   
        
    if (exchangeStrategy == 2) {   
        densitytemp *= 1 / density0mean;
        
        tempX = Dxf * densitytemp;    
        tempY = Dyf * densitytemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialLength);   
        
        modelDerivativeXtemp += tempX;
        modelDerivativeYtemp += tempY;
    }
    
    dataMisfit.setModelDerivativeX(modelDerivativeXtemp);
    dataMisfit.setModelDerivativeY(modelDerivativeYtemp);
}

/*! \brief Function for calculating the cross gradient 
 *
 \param dataMisfitEM EM dataMisfit to store the modelEM derivatives.
 \param model model parameters.
 \param derivatives derivatives
 \param config config handle
 \param modelTaper2DJoint Taper to provide transform matrix
 \param workflow workflow
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, KITGPI::Configuration::Configuration config, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityStemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempXY;
    scai::lama::DenseVector<ValueType> tempYX;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(config, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(config, "NY");
    scai::IndexType spatialLength = config.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = config.get<IndexType>("exchangeStrategy");
    
    /* Get references to required derivatives matrices */
    scai::lama::CSRSparseMatrix<ValueType> Dxf;
    scai::lama::CSRSparseMatrix<ValueType> Dyf;
    ValueType DT = config.get<ValueType>("DT");
    Dxf = derivatives.getDxf();
    Dyf = derivatives.getDyf();
    Dxf.scale(1.0 / DT);
    Dyf.scale(1.0 / DT);
                      
    scai::hmemo::ContextPtr ctx = model.getDensity().getContextPtr();
    scai::dmemo::DistributionPtr dist = model.getDensity().getDistributionPtr();    
                    
    densitytemp = model.getDensity();  
    velocityStemp = model.getVelocityS(); 
        
    // store the mean value of model parameters for weighting the gradient in summing
    if (workflow.workflowStage + workflow.iteration == 0) {
        density0mean = densitytemp.sum() / densitytemp.size();
        velocityS0mean = velocityStemp.sum() / velocityStemp.size();
    }
    
    velocityStemp *= 1 / velocityS0mean;
    
    modelDerivativeXtemp = Dxf * velocityStemp;    
    modelDerivativeYtemp = Dyf * velocityStemp;
    
    KITGPI::Common::applyMedianFilterTo2DVector(modelDerivativeXtemp, NX, NY, spatialLength);
    KITGPI::Common::applyMedianFilterTo2DVector(modelDerivativeYtemp, NX, NY, spatialLength);       
                                                   
    if (workflow.getInvertForVs() && exchangeStrategy != 0) {  
        // cross gradient of vs and EM model   
        modelTaper2DJoint.applyGradientTransform2to1(tempYX, dataMisfitEM.getModelDerivativeX());  
        modelTaper2DJoint.applyGradientTransform2to1(tempXY, dataMisfitEM.getModelDerivativeY());  
        tempYX *= modelDerivativeYtemp;   
        tempXY *= modelDerivativeXtemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);   
        
        velocityS = tempYX - tempXY;   
    } else {
        this->initParameterisation(velocityS, ctx, dist, 0.0);
    }       
    
    if (workflow.getInvertForDensity()) { 
        // cross gradient of density and vs   
        densitytemp *= 1 / density0mean;
        
        tempXY = Dxf * densitytemp;    
        tempYX = Dyf * densitytemp; 
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);          
        
        tempYX *= modelDerivativeXtemp;   
        tempXY *= modelDerivativeYtemp;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);   
        
        density = tempYX - tempXY;    
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }        
}

/*! \brief Function for calculating the derivatives of cross gradient 
 *
 \param dataMisfitEM EM dataMisfit to store the modelEM derivatives.
 \param model model parameters.
 \param derivatives derivatives
 \param config config handle
 \param modelTaper2DJoint Taper to provide transform matrix
 \param workflow workflow
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, KITGPI::Configuration::Configuration config, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> velocityStemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempXY;
    scai::lama::DenseVector<ValueType> tempYX;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(config, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(config, "NY");
    scai::IndexType spatialLength = config.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = config.get<IndexType>("exchangeStrategy");
                  
    /* Get references to required derivatives matrices */
    scai::lama::CSRSparseMatrix<ValueType> Dxf;
    scai::lama::CSRSparseMatrix<ValueType> Dyf;
    ValueType DT = config.get<ValueType>("DT");
    Dxf = derivatives.getDxf();
    Dyf = derivatives.getDyf();
    Dxf.scale(1.0 / DT);
    Dyf.scale(1.0 / DT);
        
    scai::hmemo::ContextPtr ctx = model.getDensity().getContextPtr();
    scai::dmemo::DistributionPtr dist = model.getDensity().getDistributionPtr(); 
                                    
    if (workflow.getInvertForVs() && exchangeStrategy != 0) {    
        // derivative of cross gradient of vs and EM model with respect to vs  
        modelTaper2DJoint.applyGradientTransform2to1(tempYX, dataMisfitEM.getModelDerivativeX());  
        modelTaper2DJoint.applyGradientTransform2to1(tempXY, dataMisfitEM.getModelDerivativeY());  
        tempYX *= velocityS;   
        tempXY *= velocityS;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);   
        
        tempYX = Dyf * tempYX;   
        tempXY = Dxf * tempXY;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);
        
        velocityS = tempXY - tempYX;   
    } else {
        this->initParameterisation(velocityS, ctx, dist, 0.0);
    }   
        
    velocityStemp = model.getVelocityS(); 
    
    velocityStemp *= 1 / velocityS0mean;
    
    modelDerivativeXtemp = Dxf * velocityStemp;    
    modelDerivativeYtemp = Dyf * velocityStemp;
    
    KITGPI::Common::applyMedianFilterTo2DVector(modelDerivativeXtemp, NX, NY, spatialLength);
    KITGPI::Common::applyMedianFilterTo2DVector(modelDerivativeYtemp, NX, NY, spatialLength);       
                                   
    if (workflow.getInvertForDensity()) {
        // derivative of cross gradient of density and vs with respect to density  
        tempYX = modelDerivativeXtemp * density; 
        tempXY = modelDerivativeYtemp * density;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);
        
        tempYX = Dyf * tempYX;   
        tempXY = Dxf * tempXY;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempXY, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(tempYX, NX, NY, spatialLength);
        
        density = tempXY - tempYX;  
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }         
}

/*! \brief calculate the misfit of CrossGradient
 *
 */
template <typename ValueType>
ValueType KITGPI::Gradient::SH<ValueType>::calcCrossGradientMisfit()
{
    ValueType misfitSum = 0;
    IndexType count = 0;
    if (velocityS.l2Norm() != 0) {
        misfitSum += velocityS.l2Norm() / velocityS.size();
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
void KITGPI::Gradient::SH<ValueType>::calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPriori, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{    
    scai::lama::DenseVector<ValueType> velocitySPrioritemp;
    scai::lama::DenseVector<ValueType> densityPrioritemp;
    
    velocityS = model.getVelocityS();  
    density = model.getDensity();  
    velocitySPrioritemp = modelPriori.getVelocityS();  
    densityPrioritemp = modelPriori.getDensity();  
            
    scai::hmemo::ContextPtr ctx = velocityS.getContextPtr();
    scai::dmemo::DistributionPtr dist = velocityS.getDistributionPtr();
    
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

template class KITGPI::Gradient::SH<float>;
template class KITGPI::Gradient::SH<double>;
