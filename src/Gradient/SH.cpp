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
};

/*! \brief Function for summing the gradients of all shot domains
 *
 \param commInterShot inter shot communication pointer
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot)
{
    /*reduction between shot domains.
    each shot domain may have a different distribution of (gradient) vectors. This happens if geographer is used (different result for dist on each shot domain even for homogenous architecture) or on heterogenous architecture. In this case even the number of processes on each domain can vary. Therfore it is necessary that only one process per shot domain communicates all data.
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
}

/*! \brief If stream configuration is used, get a pershot model from the big model
 \param model model
 \param gradientPerShot gradient per shot
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::sumGradientPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth)
{
    auto distBig = density.getDistributionPtr();
    auto dist = gradientPerShot.getDensity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> recoverMatrix = model.getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd));
    recoverMatrix.assignTranspose(recoverMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = model.getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd), boundaryWidth);
    scai::lama::SparseVector<ValueType> recoverVector;
    recoverVector = 1.0 - eraseVector;
    
    scai::lama::DenseVector<ValueType> temp;
    
    temp = recoverMatrix * gradientPerShot.getVelocityS(); //transform pershot into big model
    temp *= recoverVector;
    velocityS *= eraseVector;
    velocityS += temp; //take over the values

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
void KITGPI::Gradient::SH<ValueType>::scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{    
    ValueType maxValue = 1;      
    
    IndexType scaleGradient = config.get<IndexType>("scaleGradient");    
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
}

template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityStemp;
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
    
    dataMisfit.setModelDerivativeX(modelDerivativeXtemp);
    dataMisfit.setModelDerivativeY(modelDerivativeYtemp);
}

template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityStemp;
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
}

template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityStemp;
    scai::lama::DenseVector<ValueType> modelEMDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelEMDerivativeYtemp;
    scai::lama::DenseVector<ValueType> temp;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialFDorder = configEM.get<IndexType>("spatialFDorder");
        
    densitytemp = this->getDensity();
    velocityStemp = this->getVelocityS();    
        
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
}

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

/*! \brief Apply a median filter to filter the extrame value of the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::SH<ValueType>::applyMedianFilter(KITGPI::Configuration::Configuration config)
{
    scai::lama::DenseVector<ValueType> density_temp;
    scai::lama::DenseVector<ValueType> velocityS_temp;
    scai::lama::DenseVector<ValueType> porosity_temp;
    scai::lama::DenseVector<ValueType> saturation_temp;
    
    density_temp = this->getDensity();
    velocityS_temp = this->getVelocityS();
    porosity_temp = this->getPorosity();
    saturation_temp = this->getSaturation();

    scai::IndexType NZ = config.get<IndexType>("NZ");
    if (NZ == 1) {
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity_temp.size() / NY;
        scai::IndexType spatialFDorder = config.get<IndexType>("spatialFDorder");
            
        KITGPI::Common::applyMedianFilterTo2DVector(density_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(velocityS_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(porosity_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(saturation_temp, NX, NY, spatialFDorder);
        
        this->setDensity(density_temp);
        this->setVelocityS(velocityS_temp);
        this->setPorosity(porosity_temp);    
        this->setSaturation(saturation_temp);
    }
}

template class KITGPI::Gradient::SH<float>;
template class KITGPI::Gradient::SH<double>;
