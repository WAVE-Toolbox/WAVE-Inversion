#include "GradientSeismic.hpp"

using namespace scai;
using namespace KITGPI;
   
/*! \brief initialize an inner workflow */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::setInvertForParameters(std::vector<bool> setInvertForParameters)
{
    workflowInner.isSeismic = true;
    workflowInner.setInvertForParameters(setInvertForParameters);
}

/*! \brief calculate the derivative of density with respect to porosity
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientSeismic<ValueType>::getDensityDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{   
    scai::lama::DenseVector<ValueType> rho_satDePorosity;
    scai::lama::DenseVector<ValueType> saturationtemp;   
    scai::lama::DenseVector<ValueType> temp1;   
    
    saturationtemp = model.getSaturation(); 

    // derivative of density 
    // \pdv{\rho_{sat}}{\phi}=S_w \rho_{w} + \left(1-S_w\right) \rho_{a} - \rho_{ma}
    ValueType const DensityWater = model.getDensityWater();
    ValueType const DensityAir = model.getDensityAir();
    rho_satDePorosity = saturationtemp * DensityWater;  
    temp1 = 1 - saturationtemp;  
    temp1 *= DensityAir; 
    rho_satDePorosity += temp1;
    rho_satDePorosity -= model.getDensityRockMatrix();  
    
    return (rho_satDePorosity);
}

/*! \brief calculate the derivative of Mu_sat with respect to porosity
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientSeismic<ValueType>::getMu_satDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{   
    scai::lama::DenseVector<ValueType> mu_satDePorosity;
    scai::lama::DenseVector<ValueType> betaDePorosity;
    scai::lama::DenseVector<ValueType> mu_ma;
    
    mu_ma = model.getShearModulusRockMatrix();
    betaDePorosity = this->getBiotCoefficientDePorosity(model);

    // derivative of mu_sat
    // \pdv{\mu_{sat}}{\phi}=-\mu_{ma} \pdv{\beta}{\phi}
    mu_satDePorosity = -mu_ma;
    mu_satDePorosity *= betaDePorosity;   
    
    return (mu_satDePorosity);
}

/*! \brief calculate the derivative of K_sat with respect to porosity
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientSeismic<ValueType>::getK_satDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{   
    scai::lama::DenseVector<ValueType> K_satDePorosity;
    scai::lama::DenseVector<ValueType> betaDePorosity;
    scai::lama::DenseVector<ValueType> Kf;
    scai::lama::DenseVector<ValueType> K_ma;
    scai::lama::DenseVector<ValueType> M;
    scai::lama::DenseVector<ValueType> beta;
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    
    K_ma = model.getBulkModulusRockMatrix();
    beta = model.getBiotCoefficient();
    betaDePorosity = this->getBiotCoefficientDePorosity(model);
    Kf = model.getBulkModulusKf();
    M = model.getBulkModulusM();

    // derivative of K_sat
    // \pdv{K_{sat}}{\phi}=\left(-\frac{\beta^2 M^2}{K_{ma}} + 2\beta M - K_{ma}\right) \pdv{\beta}{\phi} - \beta^2 M^2 \left( \frac{1}{K_{f}} - \frac{1}{K_{ma}}\right)
    temp1 = 1 / K_ma;
    Kf = 1 / Kf;
    K_satDePorosity = Kf - temp1;
    temp1 = beta * M;       
    temp2 = 2 * temp1;
    temp1 *= temp1;
    K_satDePorosity *= temp1;
    temp1 /= K_ma;
    temp1 += K_ma;
    temp2 -= temp1;
    temp2 *= betaDePorosity;
    K_satDePorosity = temp2 - K_satDePorosity;  
    
    return (K_satDePorosity);
}

/*! \brief calculate the derivative of BiotCoefficient beta with respect to porosity
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientSeismic<ValueType>::getBiotCoefficientDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    scai::lama::DenseVector<ValueType> beta;
    scai::lama::DenseVector<ValueType> betaDePorosity;
    scai::lama::DenseVector<ValueType> mask;
    
    beta = model.getBiotCoefficient();
    
    scai::dmemo::DistributionPtr dist = beta.getDistributionPtr();    
    betaDePorosity.allocate(dist);
    ValueType const CriticalPorosity = model.getCriticalPorosity();
    betaDePorosity.setScalar(1 / CriticalPorosity);
    
    mask = 1 - beta;
    mask.unaryOp(mask, common::UnaryOp::SIGN);
    betaDePorosity *= mask;
    
    return (betaDePorosity);
}

/*! \brief calculate the derivative of density with respect to saturation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientSeismic<ValueType>::getDensityDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{   
    scai::lama::DenseVector<ValueType> densityDeSaturation;
    scai::lama::DenseVector<ValueType> porositytemp;
    ValueType tempValue; 
    
    porositytemp = model.getPorosity();  
    
    // derivative of density
    // \pdv{\rho_{sat}}{S_w}=\phi\left( \rho_{w} - \rho_{a}\right)
    ValueType const DensityWater = model.getDensityWater();
    ValueType const DensityAir = model.getDensityAir();
    tempValue = DensityWater - DensityAir;
    densityDeSaturation = porositytemp * tempValue;
    
    return (densityDeSaturation);
}

/*! \brief calculate the derivative of K_sat with respect to saturation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientSeismic<ValueType>::getK_satDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{   
    scai::lama::DenseVector<ValueType> K_satDeSaturation;
    scai::lama::DenseVector<ValueType> beta;
    scai::lama::DenseVector<ValueType> M;
    scai::lama::DenseVector<ValueType> porositytemp; 
    scai::lama::DenseVector<ValueType> temp1;
    ValueType tempValue;
    
    porositytemp = model.getPorosity();      
    beta = model.getBiotCoefficient();
    M = model.getBulkModulusM();
    
    // derivative of K_sat
    // \pdv{K_{sat}}{S_w}=-\phi \beta^2 M^2 \left(\frac{1}{K_{w}} - \frac{1}{K_{a}}\right)
    ValueType const BulkModulusWater = model.getBulkModulusWater();
    ValueType const BulkModulusAir = model.getBulkModulusAir();
    tempValue = -(1 / BulkModulusWater - 1 / BulkModulusAir);
    temp1 = beta * M;
    temp1 *=temp1;
    K_satDeSaturation = porositytemp * temp1;
    K_satDeSaturation *= tempValue;
    
    return (K_satDeSaturation);
}

/*! \brief calculate Gaussian kernel
\param model model
\param modelCoordinates coordinate class object of the model
\param FCmax max frequency
*/
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::calcGaussianKernel(scai::dmemo::CommunicatorPtr commAll, KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Configuration::Configuration config, ValueType FCmax)
{
    scai::IndexType smoothGradient = config.get<IndexType>("smoothGradient");
    if (smoothGradient == 2 || smoothGradient == 3) {
        double start_t = common::Walltime::get();
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity.size() / NY; // NX is different in stream configuration
        ValueType DH = config.get<ValueType>("DH");
        scai::lama::DenseVector<ValueType> velocity; 
        velocity = model.getVelocityS();
        ValueType velocityMean = velocity.sum() / velocity.size(); 
        
        KITGPI::Common::calcGaussianKernelFor2DVector(porosity, GaussianKernel, ksize, NX, NY, DH, velocityMean, FCmax);
        
        double end_t = common::Walltime::get();
        HOST_PRINT(commAll, "\nCalculate Gaussian kernel with matrix size = " << NX*NY << " and kernel size = " << ksize << " (" << ksize*DH << " m) in " << end_t - start_t << " sec.\n");
    }
}

/*! \brief Get const reference to density parameter
 * 
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getDensity()
{

    return (density);
}

/*! \brief Get const reference to density parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getDensity() const
{
    return (density);
}

/*! \brief Set density  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::setDensity(scai::lama::Vector<ValueType> const &setDensity)
{
    density = setDensity;
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getVelocityP()
{
    return (velocityP);
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getVelocityP() const
{
    return (velocityP);
}

/*! \brief Set P-wave velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP)
{
    velocityP = setVelocityP;
}

/*! \brief Get const reference to S-wave velocity
 * 
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getVelocityS()
{
    return (velocityS);
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getVelocityS() const
{
    return (velocityS);
}

/*! \brief Set S-wave velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS)
{
    velocityS = setVelocityS;
}

/*! \brief Get const reference to tauP
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getTauP()
{
    return (tauP);
}

/*! \brief Get const reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getTauP() const
{
    return (tauP);
}

/*! \brief Set tauP parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::setTauP(scai::lama::Vector<ValueType> const &setTauP)
{
    tauP = setTauP;
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getTauS()
{
    return (tauS);
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getTauS() const
{
    return (tauS);
}

/*! \brief Set tauS parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::setTauS(scai::lama::Vector<ValueType> const &setTauS)
{
    tauS = setTauS;
}


/* EM */
/*! \brief apply parameterisation of model parameters */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::applyParameterisation(ValueType &modelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation)
{    
    COMMON_THROWEXCEPTION("There is no applyParameterisation in the GradientSeismic.")
}

/*! \brief execute parameterisation of model parameters */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::applyParameterisation(scai::lama::DenseVector<ValueType> &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation)
{    
    COMMON_THROWEXCEPTION("There is no applyParameterisation in the GradientSeismic.")
}

/*! \brief delete parameterisation of model parameters */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::deleteParameterisation(scai::lama::DenseVector<ValueType> &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation)
{     
    COMMON_THROWEXCEPTION("There is no deleteParameterisation in the GradientSeismic.")
}

/*! \brief apply parameterisation to gradient parameters */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::gradientParameterisation(scai::lama::DenseVector<ValueType> &vecGradientParameter, scai::lama::Vector<ValueType> const &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation)
{     
    COMMON_THROWEXCEPTION("There is no gradientParameterisation in the GradientSeismic.")
}

/*! \brief calculate the derivative of conductivity with respect to porosity */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientSeismic<ValueType>::getElectricConductivityDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    COMMON_THROWEXCEPTION("There is no getElectricConductivityDePorosity in the GradientSeismic.")
    scai::lama::DenseVector<ValueType> conductivityDePorosity;
    
    return conductivityDePorosity;
}

/*! \brief calculate the derivative of dielectricPermittiviy with respect to porosity */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientSeismic<ValueType>::getDielectricPermittiviyDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    COMMON_THROWEXCEPTION("There is no getDielectricPermittiviyDePorosity in the GradientSeismic.")
    scai::lama::DenseVector<ValueType> dielectricPermittiviyDePorosity;
    
    return dielectricPermittiviyDePorosity;
}

/*! \brief calculate the derivative of conductivity with respect to saturation */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientSeismic<ValueType>::getElectricConductivityDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    COMMON_THROWEXCEPTION("There is no getElectricConductivityDeSaturation in the GradientSeismic.")
    scai::lama::DenseVector<ValueType> conductivityDeSaturation;
    
    return conductivityDeSaturation;
}

/*! \brief calculate the derivative of conductivity with respect to saturation */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientSeismic<ValueType>::getDielectricPermittiviyDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    COMMON_THROWEXCEPTION("There is no getDielectricPermittiviyDeSaturation in the GradientSeismic.")
    scai::lama::DenseVector<ValueType> dielectricPermittiviyDeSaturation;
    
    return dielectricPermittiviyDeSaturation;
}

/*! \brief Get const reference to electricConductivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getElectricConductivity()
{
    COMMON_THROWEXCEPTION("There is no getElectricConductivity in the GradientSeismic.")
    return (electricConductivity);
}

/*! \brief Get const reference to electricConductivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getElectricConductivity() const
{
    COMMON_THROWEXCEPTION("There is no getElectricConductivity in the GradientSeismic.")
    return (electricConductivity);
}

/*! \brief Set electricConductivity
 */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::setElectricConductivity(scai::lama::Vector<ValueType> const &setElectricConductivity)
{
    COMMON_THROWEXCEPTION("There is no setElectricConductivity in the GradientSeismic.")
    electricConductivity = setElectricConductivity;
}

/*! \brief Get const reference to dielectricPermittivity
 * 
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getDielectricPermittivity()
{
    COMMON_THROWEXCEPTION("There is no getDielectricPermittivity in the GradientSeismic.")
    return (dielectricPermittivity);
}

/*! \brief Get const reference to dielectricPermittivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getDielectricPermittivity() const
{
    COMMON_THROWEXCEPTION("There is no getDielectricPermittivity in the GradientSeismic.")
    return (dielectricPermittivity);
}

/*! \brief Set dielectricPermittivity
 */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::setDielectricPermittivity(scai::lama::Vector<ValueType> const &setDielectricPermittivity)
{
    COMMON_THROWEXCEPTION("There is no setDielectricPermittivity in the GradientSeismic.")
    dielectricPermittivity = setDielectricPermittivity;
}

/*! \brief Get const reference to tauElectricConductivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getTauElectricConductivity()
{
    COMMON_THROWEXCEPTION("There is no getTauElectricConductivity in the GradientSeismic.")
    return (tauElectricConductivity);
}

/*! \brief Get const reference to tauElectricConductivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getTauElectricConductivity() const
{
    COMMON_THROWEXCEPTION("There is no getTauElectricConductivity in the GradientSeismic.")
    return (tauElectricConductivity);
}

/*! \brief Set tauElectricConductivity parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::setTauElectricConductivity(scai::lama::Vector<ValueType> const &setTauElectricConductivity)
{
    COMMON_THROWEXCEPTION("There is no setTauElectricConductivity in the GradientSeismic.")
    tauElectricConductivity = setTauElectricConductivity;
}

/*! \brief Get const reference to tauDielectricPermittivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getTauDielectricPermittivity()
{
    COMMON_THROWEXCEPTION("There is no getTauDielectricPermittivity in the GradientSeismic.")
    return (tauDielectricPermittivity);
}

/*! \brief Get const reference to tauDielectricPermittivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientSeismic<ValueType>::getTauDielectricPermittivity() const
{
    COMMON_THROWEXCEPTION("There is no getTauDielectricPermittivity in the GradientSeismic.")
    return (tauDielectricPermittivity);
}

/*! \brief Set tauDielectricPermittivity parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientSeismic<ValueType>::setTauDielectricPermittivity(scai::lama::Vector<ValueType> const &setTauDielectricPermittivity)
{
    COMMON_THROWEXCEPTION("There is no setTauDielectricPermittivity in the GradientSeismic.")
    tauDielectricPermittivity = setTauDielectricPermittivity;
}

template class KITGPI::Gradient::GradientSeismic<float>;
template class KITGPI::Gradient::GradientSeismic<double>;
