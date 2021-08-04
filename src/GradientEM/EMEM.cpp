#include "EMEM.hpp"
#include <scai/dmemo/SingleDistribution.hpp>

using namespace scai;

/*! \brief Constructor that is using the Configuration class
 *
 \param ctx Context for the Calculation
 \param distEM Distribution
 */
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType>::EMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM)
{
    equationTypeEM = "emem";
    init(ctx, distEM, 0.0, 0.0, 0.0, 0.0);
}

/*! \brief Initialisation with zeros
 *
 \param ctx Context for the Calculation
 \param distEM Distribution
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM)
{
    init(ctx, distEM, 0.0, 0.0, 0.0, 0.0);
}

/*! \brief Initialisation that is generating a homogeneous modelEM
 *
 *  Generates a homogeneous modelEM, which will be initialized by the two given scalar values.
 \param ctx Context
 \param distEM Distribution
 \param electricConductivity_const electricConductivity given as Scalar
 \param dielectricPermittivity_const dielectricPermittivity given as Scalar
 \param porosity_const porosity given as Scalar
 \param saturation_const saturation given as Scalar
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, ValueType electricConductivity_const, ValueType dielectricPermittivity_const, ValueType porosity_const, ValueType saturation_const)
{
    this->initParameterisation(electricConductivity, ctx, distEM, electricConductivity_const);
    this->initParameterisation(dielectricPermittivity, ctx, distEM, dielectricPermittivity_const);
    this->initParameterisation(porosity, ctx, distEM, porosity_const);
    this->initParameterisation(saturation, ctx, distEM, saturation_const);
    this->initParameterisation(reflectivity, ctx, distEM, 0.0);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType>::EMEM(const EMEM &rhs)
{
    equationTypeEM = rhs.equationTypeEM;
    electricConductivity = rhs.electricConductivity;
    dielectricPermittivity = rhs.dielectricPermittivity;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    reflectivity = rhs.reflectivity;
}

/*! \brief Set all parameter to zero.
*/
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::resetGradient()
{
    this->resetParameter(electricConductivity);
    this->resetParameter(dielectricPermittivity);
    this->resetParameter(porosity);
    this->resetParameter(saturation);
    this->resetParameter(reflectivity);
}

/*! \brief Write modelEM to an external file
 *
 \param filename For the electricConductivity ".sigmaEM.mtx" and for dielectricPermittivity ".epsilonEMr.mtx" is added.
 \param fileFormat format of output file
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::write(std::string filename, IndexType fileFormat, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) const
{
    if (workflowEM.getInvertForSigmaEM()) {
        std::string filenameElectricConductivity = filename + ".sigmaEM";
        this->writeParameterisation(electricConductivity, filenameElectricConductivity, fileFormat);
    }
    
    if (workflowEM.getInvertForEpsilonEM()) {
        std::string filenameDielectricPermittivity = filename + ".epsilonEMr";
        this->writeParameterisation(dielectricPermittivity, filenameDielectricPermittivity, fileFormat);
    }
    
    if (workflowEM.getInvertForPorosity()) { 
        std::string filenamePorosity = filename + ".porosity";
        this->writeParameterisation(porosity, filenamePorosity, fileFormat);
    }
        
    if (workflowEM.getInvertForSaturation()) { 
        std::string filenameSaturation = filename + ".saturation";
        this->writeParameterisation(saturation, filenameSaturation, fileFormat);
    }
    
    if (workflowEM.getInvertForReflectivity()) { 
        std::string filenameReflectivity = filename + ".reflectivity";
        this->writeParameterisation(reflectivity, filenameReflectivity, fileFormat);
    }
};

/*! \brief Get equationTypeEM (emem)
 */
template <typename ValueType>
std::string KITGPI::Gradient::EMEM<ValueType>::getEquationType() const
{
    return (equationTypeEM);
}

/*! \brief Get reference to tauElectricConductivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::EMEM<ValueType>::getTauElectricConductivity()
{
    COMMON_THROWEXCEPTION("There is no tauElectricConductivity in an emem modelling")
    return (tauElectricConductivity);
}

/*! \brief Get reference to tauDielectricPermittivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::EMEM<ValueType>::getTauDielectricPermittivity()
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivity in an emem modelling")
    return (tauDielectricPermittivity);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Gradient::EMEM<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an emem modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Gradient::EMEM<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an emem modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType> KITGPI::Gradient::EMEM<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Gradient::EMEM<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief Free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType> operator*(ValueType lhs, KITGPI::Gradient::EMEM<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType> &KITGPI::Gradient::EMEM<ValueType>::operator*=(ValueType const &rhs)
{
    if (workflowInner.getInvertForSigmaEM()) 
        electricConductivity *= rhs;
    if (workflowInner.getInvertForEpsilonEM()) 
        dielectricPermittivity *= rhs;
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
KITGPI::Gradient::EMEM<ValueType> KITGPI::Gradient::EMEM<ValueType>::operator+(KITGPI::Gradient::EMEM<ValueType> const &rhs)
{
    KITGPI::Gradient::EMEM<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType> &KITGPI::Gradient::EMEM<ValueType>::operator+=(KITGPI::Gradient::EMEM<ValueType> const &rhs)
{
    electricConductivity += rhs.electricConductivity;
    dielectricPermittivity += rhs.dielectricPermittivity;
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
KITGPI::Gradient::EMEM<ValueType> KITGPI::Gradient::EMEM<ValueType>::operator-(KITGPI::Gradient::EMEM<ValueType> const &rhs)
{
    KITGPI::Gradient::EMEM<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType> &KITGPI::Gradient::EMEM<ValueType>::operator-=(KITGPI::Gradient::EMEM<ValueType> const &rhs)
{
    electricConductivity -= rhs.electricConductivity;
    dielectricPermittivity -= rhs.dielectricPermittivity;
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
KITGPI::Gradient::EMEM<ValueType> &KITGPI::Gradient::EMEM<ValueType>::operator=(KITGPI::Gradient::EMEM<ValueType> const &rhs)
{
    electricConductivity = rhs.electricConductivity;
    dielectricPermittivity = rhs.dielectricPermittivity;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    reflectivity = rhs.reflectivity;

    return *this;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::assign(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    electricConductivity = rhs.getElectricConductivity();
    dielectricPermittivity = rhs.getDielectricPermittivity();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    reflectivity = rhs.getReflectivity();
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::minusAssign(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    electricConductivity -= rhs.getElectricConductivity();
    dielectricPermittivity -= rhs.getDielectricPermittivity();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    reflectivity -= rhs.getReflectivity();
}

/*! \brief Function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::plusAssign(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    electricConductivity += rhs.getElectricConductivity();
    dielectricPermittivity += rhs.getDielectricPermittivity();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    reflectivity += rhs.getReflectivity();
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::timesAssign(ValueType const &rhs)
{
    if (workflowInner.getInvertForSigmaEM()) 
        electricConductivity *= rhs;
    if (workflowInner.getInvertForEpsilonEM()) 
        dielectricPermittivity *= rhs;
    if (workflowInner.getInvertForPorosity()) 
        porosity *= rhs;
    if (workflowInner.getInvertForSaturation()) 
        saturation *= rhs;
    if (workflowInner.getInvertForReflectivity()) 
        reflectivity *= rhs;
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::timesAssign(scai::lama::Vector<ValueType> const &rhs)
{
    electricConductivity *= rhs;
    dielectricPermittivity *= rhs;
    porosity *= rhs;
    saturation *= rhs;
    reflectivity *= rhs;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract modelEM.
 \param rhs Abstract gradientEM which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> &lhs, KITGPI::Gradient::GradientEM<ValueType> const &rhs)
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
    } else {
        scai::lama::DenseVector<ValueType> electricConductivitytemp;
        scai::lama::DenseVector<ValueType> dielectricPermittivitytemp; 
        ValueType const DielectricPermittivityVacuum = lhs.getDielectricPermittivityVacuum();
        ValueType const ElectricConductivityReference = lhs.getElectricConductivityReference();    
        
        if (workflowInner.getInvertForSigmaEM()) {
            electricConductivitytemp = lhs.getElectricConductivity();  
            this->applyParameterisation(electricConductivitytemp, ElectricConductivityReference, lhs.getParameterisation()); 
            electricConductivitytemp -= rhs.getElectricConductivity();    
            this->deleteParameterisation(electricConductivitytemp, ElectricConductivityReference, lhs.getParameterisation());  
            lhs.setElectricConductivity(electricConductivitytemp);
        }
        if (workflowInner.getInvertForEpsilonEM()) {
            dielectricPermittivitytemp = lhs.getDielectricPermittivity(); 
            this->applyParameterisation(dielectricPermittivitytemp, DielectricPermittivityVacuum, lhs.getParameterisation());  
            dielectricPermittivitytemp -= rhs.getDielectricPermittivity(); 
            this->deleteParameterisation(dielectricPermittivitytemp, DielectricPermittivityVacuum, lhs.getParameterisation());
            lhs.setDielectricPermittivity(dielectricPermittivitytemp);
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
void KITGPI::Gradient::EMEM<ValueType>::sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot)
{
    /*reduction between shot domains.
    each shot domain may have a different distribution of (gradient) vectors. This happens if geographer is used (different result for dist on each shot domain even for homogenous architecture) or on heterogenous architecture. In this case even the number of processes on each domain can vary. Therfore it is necessary that only one process per shot domain communicates all data.
    */
    
    //get information from distributed vector
    auto size = electricConductivity.size();
    auto distEM = electricConductivity.getDistributionPtr();
    auto comm = distEM->getCommunicatorPtr();
    
    // create single distribution, only master process owns the complete vector (no distribution).
    
    int shotMaster = 0;
    auto singleDist = std::make_shared<dmemo::SingleDistribution>( size, comm, shotMaster );
    
    //redistribute vector to master process
    // (this may cause memory issues for big models)
    electricConductivity.redistribute(singleDist);    
    //reduce local array (size of local array is !=0 only for master process)
    commInterShot->sumArray(electricConductivity.getLocalValues());    
    //redistribute vector to former partition
    electricConductivity.redistribute(distEM);
        
    dielectricPermittivity.redistribute(singleDist);
    commInterShot->sumArray(dielectricPermittivity.getLocalValues());
    dielectricPermittivity.redistribute(distEM);
    
    porosity.redistribute(singleDist);
    commInterShot->sumArray(porosity.getLocalValues());
    porosity.redistribute(distEM);
    
    saturation.redistribute(singleDist);
    commInterShot->sumArray(saturation.getLocalValues());
    saturation.redistribute(distEM);
    
    reflectivity.redistribute(singleDist);
    commInterShot->sumArray(reflectivity.getLocalValues());
    reflectivity.redistribute(distEM);
}

/*! \brief If stream configuration is used, set a gradient per shot into the big gradient
 \param model model
 \param gradientPerShot gradient per shot
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::sumGradientPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelEM, KITGPI::Gradient::GradientEM<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth)
{
    auto distBig = dielectricPermittivity.getDistributionPtr();
    auto dist = gradientPerShot.getDielectricPermittivity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> recoverMatrix = modelEM.getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd));
    recoverMatrix.assignTranspose(recoverMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = modelEM.getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd), boundaryWidth);
    scai::lama::SparseVector<ValueType> recoverVector;
    recoverVector = 1.0 - eraseVector;
    
    scai::lama::DenseVector<ValueType> temp;
    
    temp = recoverMatrix * gradientPerShot.getDielectricPermittivity(); //transform pershot into big model
    temp *= recoverVector;
    dielectricPermittivity *= eraseVector;
    dielectricPermittivity += temp; //take over the values
  
    temp = recoverMatrix * gradientPerShot.getElectricConductivity(); //transform pershot into big model
    temp *= recoverVector;
    electricConductivity *= eraseVector;
    electricConductivity += temp; //take over the values
    
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

/*! \brief Function for scaling the gradients with the modelEM parameter 
 *
 \param modelEM Abstract modelEM.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::scale(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM, KITGPI::Configuration::Configuration configEM)
{    
    ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
    ValueType const ElectricConductivityReference = modelEM.getElectricConductivityReference();
    ValueType maxValue;      
    
    IndexType scaleGradient = configEM.get<IndexType>("scaleGradient");
    if (scaleGradient != 0) {
        if (workflowEM.getInvertForSigmaEM() && electricConductivity.maxNorm() != 0) {  
            if (scaleGradient == 1) {
                maxValue = modelEM.getElectricConductivity().maxNorm();
            } else if (scaleGradient == 2) {
                maxValue = configEM.get<ValueType>("upperSigmaEMTh") - configEM.get<ValueType>("lowerSigmaEMTh");
            }
            this->applyParameterisation(maxValue, ElectricConductivityReference, modelEM.getParameterisation());
            electricConductivity *= 1 / electricConductivity.maxNorm() * maxValue;
        }  
        
        if (workflowEM.getInvertForEpsilonEM() && dielectricPermittivity.maxNorm() != 0) {
            if (scaleGradient == 1) {
                maxValue = modelEM.getDielectricPermittivity().maxNorm();
            } else if (scaleGradient == 2) {
                maxValue = (configEM.get<ValueType>("upperEpsilonEMrTh") - configEM.get<ValueType>("lowerEpsilonEMrTh")) * DielectricPermittivityVacuum;
            }      
            this->applyParameterisation(maxValue, DielectricPermittivityVacuum, modelEM.getParameterisation()); 
            dielectricPermittivity *= 1 / dielectricPermittivity.maxNorm() * maxValue;
        }
        
        if (workflowEM.getInvertForPorosity() && porosity.maxNorm() != 0) {
            if (modelEM.getParameterisation() == 1 || modelEM.getParameterisation() == 2) {
                if (scaleGradient == 1) {
                    maxValue = modelEM.getPorosity().maxNorm();
                } else if (scaleGradient == 2) {
                    maxValue = configEM.getAndCatch("upperPorosityTh", 1.0) - configEM.getAndCatch("lowerPorosityTh", 0.0);
                }
            }        
            porosity *= 1 / porosity.maxNorm() * maxValue;
        }    
    
        if (workflowEM.getInvertForSaturation() && saturation.maxNorm() != 0) { 
            if (modelEM.getParameterisation() == 1 || modelEM.getParameterisation() == 2) {
                if (scaleGradient == 1) {
                    maxValue = modelEM.getSaturation().maxNorm();
                } else if (scaleGradient == 2) {
                    maxValue = configEM.getAndCatch("upperSaturationTh", 1.0) - configEM.getAndCatch("lowerSaturationTh", 0.0);
                }
            }              
            saturation *= 1 / saturation.maxNorm() * maxValue;        
        }
        
        if (workflowEM.getInvertForReflectivity() && reflectivity.maxNorm() != 0) { 
            if (scaleGradient == 1) {
                maxValue = modelEM.getReflectivity().maxNorm();
            } else if (scaleGradient == 2) {
                maxValue = configEM.getAndCatch("upperReflectivityEMrTh", 1.0) - configEM.getAndCatch("lowerReflectivityEMrTh", -1.0);
            }      
            reflectivity *= 1 / reflectivity.maxNorm() * maxValue;      
        }
    }
}

/*! \brief Function for normalizing the gradientEM
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::normalize()
{
    if (this->getNormalizeGradient()) {
        ValueType gradientMax = electricConductivity.maxNorm();
        if (gradientMax != 0)
            electricConductivity *= 1 / gradientMax;
        gradientMax = dielectricPermittivity.maxNorm();
        if (gradientMax != 0)
            dielectricPermittivity *= 1 / gradientMax;
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

/*! \brief Function for calculating the emem gradients from the cross correlation and the modelEM parameter 
 *
 \param modelEM Abstract modelEM.
 \param correlatedWavefields Abstract xCorr.
 \param DT Temporal discretization 
 \param workflowEM 
 *
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, ValueType DT, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{    
    scai::lama::DenseVector<ValueType> gradEpsilonEM;
    scai::lama::DenseVector<ValueType> temp;     
    
    gradEpsilonEM = DT * correlatedWavefields.getXcorrEpsilonEM();
    
    scai::hmemo::ContextPtr ctx = gradEpsilonEM.getContextPtr();
    scai::dmemo::DistributionPtr distEM = gradEpsilonEM.getDistributionPtr();
    
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) { 
        electricConductivity = DT * correlatedWavefields.getXcorrSigmaEM(); 
        
        temp = modelEM.getElectricConductivity();
        ValueType const ElectricConductivityReference = modelEM.getElectricConductivityReference();
        this->gradientParameterisation(electricConductivity, temp, ElectricConductivityReference, modelEM.getParameterisation());
    } else {
        this->initParameterisation(electricConductivity, ctx, distEM, 0.0);
    }
    
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation() || workflowEM.getInvertForReflectivity()) {
        dielectricPermittivity = gradEpsilonEM;
        
        temp = modelEM.getDielectricPermittivity();
        ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
        this->gradientParameterisation(dielectricPermittivity, temp, DielectricPermittivityVacuum, modelEM.getParameterisation());
    } else {
        this->initParameterisation(dielectricPermittivity, ctx, distEM, 0.0);
    }
        
    if (workflowEM.getInvertForPorosity()) {
        // porosity and water saturation        
        scai::lama::DenseVector<ValueType> conductivityDePorosity;
        scai::lama::DenseVector<ValueType> dielectricPermittiviyDePorosity;
     
        // Based on complex refractive index model (CRIM)    
        dielectricPermittiviyDePorosity = this->getDielectricPermittiviyDePorosity(modelEM);             
        
        porosity = dielectricPermittiviyDePorosity * dielectricPermittivity;
        if (modelEM.getParameterisation() == 2) {
            // Based on Archie equation
            conductivityDePorosity = this->getElectricConductivityDePorosity(modelEM);    
            conductivityDePorosity *= electricConductivity;
            porosity += conductivityDePorosity; 
        } 
    } else {
        this->initParameterisation(porosity, ctx, distEM, 0.0);
    }
    
    if (workflowEM.getInvertForSaturation()) {
        // porosity and water saturation        
        scai::lama::DenseVector<ValueType> conductivityDeSaturation;
        scai::lama::DenseVector<ValueType> dielectricPermittiviyDeSaturation;
               
        // Based on complex refractive index model (CRIM)            
        dielectricPermittiviyDeSaturation = this->getDielectricPermittiviyDeSaturation(modelEM);             
        
        saturation = dielectricPermittiviyDeSaturation * dielectricPermittivity;
        if (modelEM.getParameterisation() == 2) {
            // Based on Archie equation
            conductivityDeSaturation = this->getElectricConductivityDeSaturation(modelEM); 
            conductivityDeSaturation *= electricConductivity;
            saturation += conductivityDeSaturation; 
        }           
    } else {
        this->initParameterisation(saturation, ctx, distEM, 0.0);
    }  
    
    if (workflowEM.getInvertForReflectivity()) {
        reflectivity = -dielectricPermittivity;
    } else {
        this->initParameterisation(reflectivity, ctx, distEM, 0.0);
    }  
}

/*! \brief calculate the stabilizing functional of each model parameter */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::calcStabilizingFunctionalGradient(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelPrioriEM, KITGPI::Configuration::Configuration configEM, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{    
    scai::lama::DenseVector<ValueType> electricConductivityPrioritemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivityPrioritemp;
    ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
    ValueType const ElectricConductivityReference = modelEM.getElectricConductivityReference();
    
    electricConductivity = modelEM.getElectricConductivity();  
    dielectricPermittivity = modelEM.getDielectricPermittivity();     
    electricConductivityPrioritemp = modelPrioriEM.getElectricConductivity();  
    dielectricPermittivityPrioritemp = modelPrioriEM.getDielectricPermittivity();

    this->applyParameterisation(electricConductivity, ElectricConductivityReference, modelEM.getParameterisation());       
    this->applyParameterisation(dielectricPermittivity, DielectricPermittivityVacuum, modelEM.getParameterisation()); 
    this->applyParameterisation(electricConductivityPrioritemp, ElectricConductivityReference, modelEM.getParameterisation());       
    this->applyParameterisation(dielectricPermittivityPrioritemp, DielectricPermittivityVacuum, modelEM.getParameterisation()); 
        
    scai::hmemo::ContextPtr ctx = dielectricPermittivity.getContextPtr();
    scai::dmemo::DistributionPtr distEM = dielectricPermittivity.getDistributionPtr();
    
    if (workflowEM.getInvertForEpsilonEM()) {
        dielectricPermittivityPrioritemp = dielectricPermittivity - dielectricPermittivityPrioritemp;
        dielectricPermittivity = this->calcStabilizingFunctionalGradientPerModel(dielectricPermittivityPrioritemp, configEM, dataMisfitEM);
    } else {
        this->initParameterisation(dielectricPermittivity, ctx, distEM, 0.0);
    }
    
    if (workflowEM.getInvertForSigmaEM()) { 
        electricConductivityPrioritemp = electricConductivity - electricConductivityPrioritemp;
        electricConductivity = this->calcStabilizingFunctionalGradientPerModel(electricConductivityPrioritemp, configEM, dataMisfitEM);
    } else {
        this->initParameterisation(electricConductivity, ctx, distEM, 0.0);
    }    
}

/*! \brief Apply a median filter to filter the extrame value of the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::applyMedianFilter(KITGPI::Configuration::Configuration configEM)
{
    scai::lama::DenseVector<ValueType> sigmaEM_temp;
    scai::lama::DenseVector<ValueType> epsilonEM_temp;
    scai::lama::DenseVector<ValueType> porosity_temp;
    scai::lama::DenseVector<ValueType> saturation_temp;
    scai::lama::DenseVector<ValueType> reflectivity_temp;
    
    sigmaEM_temp = this->getElectricConductivity();
    epsilonEM_temp = this->getDielectricPermittivity();
    porosity_temp = this->getPorosity();
    saturation_temp = this->getSaturation();
    reflectivity_temp = this->getReflectivity();
    
    scai::IndexType NZ = configEM.get<IndexType>("NZ");
    if (NZ == 1) {
        scai::IndexType NY = configEM.get<IndexType>("NY");
        scai::IndexType NX = porosity_temp.size() / NY;
        scai::IndexType spatialFDorder = configEM.get<IndexType>("spatialFDorder");
        
        KITGPI::Common::applyMedianFilterTo2DVector(sigmaEM_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(epsilonEM_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(porosity_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(saturation_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(reflectivity_temp, NX, NY, spatialFDorder);
        
        this->setElectricConductivity(sigmaEM_temp);    
        this->setDielectricPermittivity(epsilonEM_temp);
        this->setPorosity(porosity_temp);    
        this->setSaturation(saturation_temp);
        this->setReflectivity(reflectivity_temp);
    }
}

template class KITGPI::Gradient::EMEM<float>;
template class KITGPI::Gradient::EMEM<double>;
