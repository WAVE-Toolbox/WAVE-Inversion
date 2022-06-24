#include "GradientEM.hpp"

using namespace scai;
using namespace KITGPI;
    
/*! \brief initialize an inner workflow */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setInvertForParameters(std::vector<bool> setInvertForParameters)
{
    workflowInner.isSeismic = false;
    workflowInner.setInvertForParameters(setInvertForParameters);
}

/*! \brief apply parameterisation of model parameters */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::applyParameterisation(ValueType &modelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation)
{    
    switch (parameterisation) {
        case 3:            // sqrt  
            // \sqrt{x/x_0}      
            modelParameter /= modelParameterReference;                 
            modelParameter = sqrt(modelParameter); 
            break;
        case 4:            // logarithmic  
            // \ln(x/x_0+1)                   
            modelParameter /= modelParameterReference;
            modelParameter += 1; // in case that modelParameter = 0
            modelParameter = log(modelParameter); 
            if (std::isnan(modelParameter) || modelParameter == std::numeric_limits<ValueType>::infinity() || -modelParameter == std::numeric_limits<ValueType>::infinity()) {
                modelParameter = 0.0;
            }
            break;
        default:          // case 0: no parameterisation 
            // x/x_0   
            modelParameter /= modelParameterReference;  
    }
}

/*! \brief execute parameterisation of model parameters */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::applyParameterisation(scai::lama::DenseVector<ValueType> &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation)
{    
    switch (parameterisation) {
        case 3:            // sqrt  
            // \sqrt{x/x_0}      
            vecModelParameter /= modelParameterReference;                 
            vecModelParameter = scai::lama::sqrt(vecModelParameter); 
            break;
        case 4:            // logarithmic  
            // \ln(x/x_0+1)                   
            vecModelParameter /= modelParameterReference;
            vecModelParameter += 1; // in case that vecModelParameter = 0
            vecModelParameter = scai::lama::log(vecModelParameter); 
            Common::replaceInvalid<ValueType>(vecModelParameter, 0.0);
            break;
        default:          // case 0: no parameterisation 
            // x/x_0   
            vecModelParameter /= modelParameterReference;  
    }
}

/*! \brief delete parameterisation of model parameters */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::deleteParameterisation(scai::lama::DenseVector<ValueType> &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation)
{     
    switch (parameterisation) {
        case 3:            // sqrt   
            vecModelParameter = scai::lama::pow(vecModelParameter, 2.0);
            vecModelParameter *= modelParameterReference; 
            break;
        case 4:            // logarithmic
            vecModelParameter = scai::lama::exp(vecModelParameter);
            vecModelParameter -= 1;  // in case that vecModelParameter = 0
            vecModelParameter *= modelParameterReference;   
            break;
        default:          // case 0: no parameterisation  
            vecModelParameter *= modelParameterReference;    
    }
}

/*! \brief apply parameterisation to gradient parameters */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::gradientParameterisation(scai::lama::DenseVector<ValueType> &vecGradientParameter, scai::lama::Vector<ValueType> const &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation)
{     
    scai::lama::DenseVector<ValueType> temp;
    switch (parameterisation) {
        case 3:            // sqrt   
            // 2 \sqrt{x x_0} \nabla f_1(x)           
            temp = vecModelParameter * modelParameterReference;
            temp = scai::lama::sqrt(temp);
            vecGradientParameter *= temp;
            vecGradientParameter *= 2; 
            break;
        case 4:            // logarithmic 
            // (x+x_0)\nabla f_1(x)
            temp = vecModelParameter + modelParameterReference;  // in case that vecModelParameter = 0
            vecGradientParameter *= temp; 
            break;
        default:          // case 0: no parameterisation    
            // x_0 \nabla f_1(x)           
            vecGradientParameter *= modelParameterReference;    
    }
}

/*! \brief calculate the derivative of conductivity with respect to porosity */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getElectricConductivityDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    scai::lama::DenseVector<ValueType> conductivityDePorosity;
    scai::lama::DenseVector<ValueType> temp;  

    // \pdv{\sigma_e}{\phi} = \frac{m}{a} \sigma_{ew} \phi^{m-1} S_w^{n}
    ValueType aArchietemp;
    ValueType mArchietemp;
    ValueType nArchietemp;
    aArchietemp = model.getArchie_a();
    mArchietemp = model.getArchie_m();
    nArchietemp = model.getArchie_n();
    // Based on Archie equation
    conductivityDePorosity = model.getElectricConductivityWater();   
    conductivityDePorosity *= mArchietemp / aArchietemp;
    temp = scai::lama::pow(model.getPorosity(), mArchietemp - 1); 
    Common::replaceInvalid<ValueType>(temp, 0.0);
    conductivityDePorosity *= temp;
    temp = scai::lama::pow(model.getSaturation(), nArchietemp); 
    conductivityDePorosity *= temp; 
    
    return conductivityDePorosity;
}

/*! \brief calculate the derivative of dielectricPermittiviy with respect to porosity */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getDielectricPermittiviyDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    scai::lama::DenseVector<ValueType> dielectricPermittiviyDePorosity;
    scai::lama::DenseVector<ValueType> temp; 
    ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum(); 
    ValueType const RelativeDielectricPermittivityWater = model.getRelativeDielectricPermittivityWater(); 
    ValueType const RelativeDielectricPermittivityVacuum = model.getRelativeDielectricPermittivityVacuum(); 
    
    // \pdv{\varepsilon_{e}}{\phi} = 2 \left[ - \sqrt{ \varepsilon_{emar} } + S_w \sqrt{ \varepsilon_{ewr} } + \left(1-S_w\right) \sqrt{ \varepsilon_{ear}} \right] \sqrt{ \varepsilon_{e} \varepsilon_{e0}}
    // Based on complex refractive index model (CRIM)    
    dielectricPermittiviyDePorosity = scai::lama::sqrt(model.getRelativeDieletricPeimittivityRockMatrix()); 
    temp = model.getSaturation();  
    temp *= sqrt(RelativeDielectricPermittivityWater);
    dielectricPermittiviyDePorosity = temp - dielectricPermittiviyDePorosity;
    temp = 1 - model.getSaturation();  
    temp *= sqrt(RelativeDielectricPermittivityVacuum);
    dielectricPermittiviyDePorosity += temp;    
    temp = model.getDielectricPermittivity();   
    temp *= DielectricPermittivityVacuum;
    temp = scai::lama::sqrt(temp); 
    dielectricPermittiviyDePorosity *= temp;   
    dielectricPermittiviyDePorosity *= 2; 
    
    return dielectricPermittiviyDePorosity;
}

/*! \brief calculate the derivative of conductivity with respect to saturation */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getElectricConductivityDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    scai::lama::DenseVector<ValueType> conductivityDeSaturation;
    scai::lama::DenseVector<ValueType> temp;  
    
    // \pdv{\sigma_e}{S_w} = \frac{n}{a} \sigma_{ew} \phi^m S_w^{n-1}
    ValueType aArchietemp;
    ValueType mArchietemp;
    ValueType nArchietemp;
    aArchietemp = model.getArchie_a();
    mArchietemp = model.getArchie_m();
    nArchietemp = model.getArchie_n();
    // Based on Archie equation
    conductivityDeSaturation = model.getElectricConductivityWater();   
    conductivityDeSaturation *= nArchietemp / aArchietemp;
    temp = scai::lama::pow(model.getPorosity(), mArchietemp); 
    conductivityDeSaturation *= temp;
    temp = scai::lama::pow(model.getSaturation(), nArchietemp - 1); 
    Common::replaceInvalid<ValueType>(temp, 0.0);
    conductivityDeSaturation *= temp; 
    
    return conductivityDeSaturation;
}

/*! \brief calculate the derivative of conductivity with respect to saturation */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getDielectricPermittiviyDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    scai::lama::DenseVector<ValueType> dielectricPermittiviyDeSaturation;
    scai::lama::DenseVector<ValueType> temp;  
    ValueType tempValue;
    ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum();
    ValueType const RelativeDielectricPermittivityWater = model.getRelativeDielectricPermittivityWater(); 
    ValueType const RelativeDielectricPermittivityVacuum = model.getRelativeDielectricPermittivityVacuum(); 
        
    // \pdv{\varepsilon_{e}}{S_w} = 2  \phi \left( \sqrt{ \varepsilon_{ewr} } - \sqrt{ \varepsilon_{ear}} \right) \sqrt{ \varepsilon_{e} \varepsilon_{e0}}
    // Based on complex refractive index model (CRIM)    
    temp = model.getDielectricPermittivity();   
    temp *= DielectricPermittivityVacuum;
    dielectricPermittiviyDeSaturation = scai::lama::sqrt(temp);         
    temp = model.getPorosity();  
    tempValue = sqrt(RelativeDielectricPermittivityWater) - sqrt(RelativeDielectricPermittivityVacuum);
    temp *= tempValue;
    dielectricPermittiviyDeSaturation *= temp;   
    dielectricPermittiviyDeSaturation *= 2;  
    
    return dielectricPermittiviyDeSaturation;
}

/*! \brief calculate Gaussian kernel
\param model model
\param modelCoordinates coordinate class object of the model
\param FCmain main frequency
*/
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::calcGaussianKernel(scai::dmemo::CommunicatorPtr commAll, KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Configuration::Configuration config)
{
    scai::IndexType smoothGradient = config.getAndCatch("smoothGradient", 0);
    if (smoothGradient != 0) {
        double start_t = common::Walltime::get();
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity.size() / NY; // NX is different in stream configuration
        ValueType DH = config.get<ValueType>("DH");
        scai::lama::DenseVector<ValueType> velocity; 
        velocity = model.getVelocityEM();
        ValueType velocityMean = velocity.sum() / velocity.size(); 
        
        KITGPI::Common::calcGaussianKernelFor2DVector(porosity, GaussianKernel, PX, PY, NX, NY, DH, velocityMean, config.get<ValueType>("CenterFrequencyCPML"), smoothGradient);
        
        double end_t = common::Walltime::get();
        HOST_PRINT(commAll, "\nCalculate Gaussian kernel with matrix size = " << (NX+PX-1)*(NY+PY-1) << " x " << (NX+PX-1)*(NY+PY-1) << " and kernel size = " << PX << " x " << PY << " (" << PX*DH << " m x " << PY*DH << " m) in " << end_t - start_t << " sec.\n");
    }
}

/*! \brief Get const reference to electricConductivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getElectricConductivity()
{
    return (electricConductivity);
}

/*! \brief Get const reference to electricConductivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getElectricConductivity() const
{
    return (electricConductivity);
}

/*! \brief Set electricConductivity
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setElectricConductivity(scai::lama::Vector<ValueType> const &setElectricConductivity)
{
    electricConductivity = setElectricConductivity;
}

/*! \brief Get const reference to dielectricPermittivity
 * 
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getDielectricPermittivity()
{
    return (dielectricPermittivity);
}

/*! \brief Get const reference to dielectricPermittivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getDielectricPermittivity() const
{
    return (dielectricPermittivity);
}

/*! \brief Set dielectricPermittivity
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setDielectricPermittivity(scai::lama::Vector<ValueType> const &setDielectricPermittivity)
{
    dielectricPermittivity = setDielectricPermittivity;
}

/*! \brief Get const reference to tauElectricConductivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauElectricConductivity()
{
    return (tauElectricConductivity);
}

/*! \brief Get const reference to tauElectricConductivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauElectricConductivity() const
{
    return (tauElectricConductivity);
}

/*! \brief Set tauElectricConductivity parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setTauElectricConductivity(scai::lama::Vector<ValueType> const &setTauElectricConductivity)
{
    tauElectricConductivity = setTauElectricConductivity;
}

/*! \brief Get const reference to tauDielectricPermittivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauDielectricPermittivity()
{
    return (tauDielectricPermittivity);
}

/*! \brief Get const reference to tauDielectricPermittivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauDielectricPermittivity() const
{
    return (tauDielectricPermittivity);
}

/*! \brief Set tauDielectricPermittivity parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setTauDielectricPermittivity(scai::lama::Vector<ValueType> const &setTauDielectricPermittivity)
{
    tauDielectricPermittivity = setTauDielectricPermittivity;
}



/* Seismic */   
/*! \brief calculate the derivative of density with respect to porosity
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientEM<ValueType>::getDensityDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{   
    COMMON_THROWEXCEPTION("There is no getDensityDePorosity in the GradientEM.")
    scai::lama::DenseVector<ValueType> rho_satDePorosity; 
    
    return (rho_satDePorosity);
}

/*! \brief calculate the derivative of Mu_sat with respect to porosity
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientEM<ValueType>::getMu_satDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{   
    COMMON_THROWEXCEPTION("There is no getMu_satDePorosity in the GradientEM.")
    scai::lama::DenseVector<ValueType> mu_satDePorosity; 
    
    return (mu_satDePorosity);
}

/*! \brief calculate the derivative of K_sat with respect to porosity
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientEM<ValueType>::getK_satDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{   
    COMMON_THROWEXCEPTION("There is no getK_satDePorosity in the GradientEM.")
    scai::lama::DenseVector<ValueType> K_satDePorosity; 
    
    return (K_satDePorosity);
}

/*! \brief calculate the derivative of BiotCoefficient beta with respect to porosity
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientEM<ValueType>::getBiotCoefficientDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    COMMON_THROWEXCEPTION("There is no getBiotCoefficientDePorosity in the GradientEM.")
    scai::lama::DenseVector<ValueType> betaDePorosity;
    
    return (betaDePorosity);
}

/*! \brief calculate the derivative of density with respect to saturation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientEM<ValueType>::getDensityDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{   
    COMMON_THROWEXCEPTION("There is no getDensityDeSaturation in the GradientEM.")
    scai::lama::DenseVector<ValueType> densityDeSaturation;
    
    return (densityDeSaturation);
}

/*! \brief calculate the derivative of K_sat with respect to saturation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientEM<ValueType>::getK_satDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{   
    COMMON_THROWEXCEPTION("There is no getK_satDeSaturation in the GradientEM.")
    scai::lama::DenseVector<ValueType> K_satDeSaturation;
    
    return (K_satDeSaturation);
}

/*! \brief Get const reference to density parameter
 * 
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getDensity()
{
    COMMON_THROWEXCEPTION("There is no getDensity in the GradientEM.")
    return (density);
}

/*! \brief Get const reference to density parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getDensity() const
{
    COMMON_THROWEXCEPTION("There is no getDensity in the GradientEM.")
    return (density);
}

/*! \brief Set density  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setDensity(scai::lama::Vector<ValueType> const &setDensity)
{
    COMMON_THROWEXCEPTION("There is no setDensity in the GradientEM.")
    density = setDensity;
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getVelocityP()
{
    COMMON_THROWEXCEPTION("There is no getVelocityP in the GradientEM.")
    return (velocityP);
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getVelocityP() const
{
    COMMON_THROWEXCEPTION("There is no getVelocityP in the GradientEM.")
    return (velocityP);
}

/*! \brief Set P-wave velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP)
{
    COMMON_THROWEXCEPTION("There is no setVelocityP in the GradientEM.")
    velocityP = setVelocityP;
}

/*! \brief Get const reference to S-wave velocity
 * 
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getVelocityS()
{
    COMMON_THROWEXCEPTION("There is no getVelocityS in the GradientEM.")
    return (velocityS);
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getVelocityS() const
{
    COMMON_THROWEXCEPTION("There is no getVelocityS in the GradientEM.")
    return (velocityS);
}

/*! \brief Set S-wave velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS)
{
    COMMON_THROWEXCEPTION("There is no setVelocityS in the GradientEM.")
    velocityS = setVelocityS;
}

/*! \brief Get const reference to tauP
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauP()
{
    COMMON_THROWEXCEPTION("There is no getTauP in the GradientEM.")
    return (tauP);
}

/*! \brief Get const reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauP() const
{
    COMMON_THROWEXCEPTION("There is no getTauP in the GradientEM.")
    return (tauP);
}

/*! \brief Set tauP parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setTauP(scai::lama::Vector<ValueType> const &setTauP)
{
    COMMON_THROWEXCEPTION("There is no setTauP in the GradientEM.")
    tauP = setTauP;
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauS()
{
    COMMON_THROWEXCEPTION("There is no getTauS in the GradientEM.")
    return (tauS);
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauS() const
{
    COMMON_THROWEXCEPTION("There is no getTauS in the GradientEM.")
    return (tauS);
}

/*! \brief Set tauS parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setTauS(scai::lama::Vector<ValueType> const &setTauS)
{
    COMMON_THROWEXCEPTION("There is no setTauS in the GradientEM.")
    tauS = setTauS;
}

template class KITGPI::Gradient::GradientEM<float>;
template class KITGPI::Gradient::GradientEM<double>;
