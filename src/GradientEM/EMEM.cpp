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
 \param conductivityEM_const conductivityEM given as Scalar
 \param dielectricPermittivityEM_const dielectricPermittivityEM given as Scalar
 \param porosity_const porosity given as Scalar
 \param saturation_const saturation given as Scalar
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, ValueType conductivityEM_const, ValueType dielectricPermittivityEM_const, ValueType porosity_const, ValueType saturation_const)
{
    this->initParameterisation(conductivityEM, ctx, distEM, conductivityEM_const);
    this->initParameterisation(dielectricPermittivityEM, ctx, distEM, dielectricPermittivityEM_const);
    this->initParameterisation(porosity, ctx, distEM, porosity_const);
    this->initParameterisation(saturation, ctx, distEM, saturation_const);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType>::EMEM(const EMEM &rhs)
{
    equationTypeEM = rhs.equationTypeEM;
    conductivityEM = rhs.conductivityEM;
    dielectricPermittivityEM = rhs.dielectricPermittivityEM;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
}

/*! \brief Set all parameter to zero.
*/
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::resetGradient()
{
    this->resetParameter(conductivityEM);
    this->resetParameter(dielectricPermittivityEM);
    this->resetParameter(porosity);
    this->resetParameter(saturation);
}

/*! \brief Write modelEM to an external file
 *
 \param filename For the conductivityEM ".sigmaEM.mtx" and for dielectricPermittivityEM ".epsilonEMr.mtx" is added.
 \param fileFormat format of output file
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::write(std::string filename, IndexType fileFormat, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) const
{
    if (workflowEM.getInvertForSigmaEM()) {
        std::string filenameConductivityEM = filename + ".sigmaEM";
        this->writeParameterisation(conductivityEM, filenameConductivityEM, fileFormat);
    }
    
    if (workflowEM.getInvertForEpsilonEM()) {
        std::string filenameEpsilonEM = filename + ".epsilonEMr";
        this->writeParameterisation(dielectricPermittivityEM, filenameEpsilonEM, fileFormat);
    }
    
    if (workflowEM.getInvertForPorosity()) { 
        std::string filenamePorosity = filename + ".porosity";
        this->writeParameterisation(porosity, filenamePorosity, fileFormat);
    }
        
    if (workflowEM.getInvertForSaturation()) { 
        std::string filenameSaturation = filename + ".saturation";
        this->writeParameterisation(saturation, filenameSaturation, fileFormat);
    }
};

/*! \brief Get equationTypeEM (emem)
 */
template <typename ValueType>
std::string KITGPI::Gradient::EMEM<ValueType>::getEquationType() const
{
    return (equationTypeEM);
}

/*! \brief Get reference to tauConductivityEM
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::EMEM<ValueType>::getTauConductivityEM()
{
    COMMON_THROWEXCEPTION("There is no tauConductivityEM in an emem modelling")
    return (tauConductivityEM);
}

/*! \brief Get reference to tauDielectricPermittivityEM
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::EMEM<ValueType>::getTauDielectricPermittivityEM()
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivityEM in an emem modelling")
    return (tauDielectricPermittivityEM);
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
        conductivityEM *= rhs;
    if (workflowInner.getInvertForEpsilonEM()) 
        dielectricPermittivityEM *= rhs;
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
    conductivityEM += rhs.conductivityEM;
    dielectricPermittivityEM += rhs.dielectricPermittivityEM;
    porosity += rhs.porosity;
    saturation += rhs.saturation;

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
    conductivityEM -= rhs.conductivityEM;
    dielectricPermittivityEM -= rhs.dielectricPermittivityEM;
    porosity -= rhs.porosity;
    saturation -= rhs.saturation;

    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Gradient which is copied.
 */
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType> &KITGPI::Gradient::EMEM<ValueType>::operator=(KITGPI::Gradient::EMEM<ValueType> const &rhs)
{
    conductivityEM = rhs.conductivityEM;
    dielectricPermittivityEM = rhs.dielectricPermittivityEM;
    porosity = rhs.porosity;
    saturation = rhs.saturation;

    return *this;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::assign(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    conductivityEM = rhs.getConductivityEM();
    dielectricPermittivityEM = rhs.getDielectricPermittivityEM();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::minusAssign(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    conductivityEM -= rhs.getConductivityEM();
    dielectricPermittivityEM -= rhs.getDielectricPermittivityEM();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
}

/*! \brief Function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::plusAssign(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    conductivityEM += rhs.getConductivityEM();
    dielectricPermittivityEM += rhs.getDielectricPermittivityEM();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::timesAssign(ValueType const &rhs)
{
    if (workflowInner.getInvertForSigmaEM()) 
        conductivityEM *= rhs;
    if (workflowInner.getInvertForEpsilonEM()) 
        dielectricPermittivityEM *= rhs;
    if (workflowInner.getInvertForPorosity()) 
        porosity *= rhs;
    if (workflowInner.getInvertForSaturation()) 
        saturation *= rhs;
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::timesAssign(scai::lama::Vector<ValueType> const &rhs)
{
    conductivityEM *= rhs;
    dielectricPermittivityEM *= rhs;
    porosity *= rhs;
    saturation *= rhs;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract modelEM.
 \param rhs Abstract gradientEM which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> &lhs, KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{        
    if (lhs.getParameterisation() == 1 || lhs.getParameterisation() == 2) { 
        scai::lama::DenseVector<ValueType> temp;  
        if (workflowInner.getInvertForPorosity()) {
            temp = lhs.getPorosity() - rhs.getPorosity();
            lhs.setPorosity(temp);
        }
        if (workflowInner.getInvertForSaturation()) {
            temp = lhs.getSaturation() - rhs.getSaturation();
            lhs.setSaturation(temp);     
        }
    } else {
        scai::lama::DenseVector<ValueType> conductivityEMtemp;
        scai::lama::DenseVector<ValueType> dielectricPermittivityEMtemp; 
        ValueType const DielectricPermittivityVacuum = lhs.getDielectricPermittivityVacuum();
        ValueType const ConductivityReference = lhs.getConductivityReference();    
        
        if (workflowInner.getInvertForSigmaEM()) {
            conductivityEMtemp = lhs.getConductivityEM();  
            this->applyParameterisation(conductivityEMtemp, ConductivityReference, lhs.getParameterisation()); 
            conductivityEMtemp -= rhs.getConductivityEM();    
            this->deleteParameterisation(conductivityEMtemp, ConductivityReference, lhs.getParameterisation());  
            lhs.setConductivityEM(conductivityEMtemp);
        }
        if (workflowInner.getInvertForEpsilonEM()) {
            dielectricPermittivityEMtemp = lhs.getDielectricPermittivityEM(); 
            this->applyParameterisation(dielectricPermittivityEMtemp, DielectricPermittivityVacuum, lhs.getParameterisation());  
            dielectricPermittivityEMtemp -= rhs.getDielectricPermittivityEM(); 
            this->deleteParameterisation(dielectricPermittivityEMtemp, DielectricPermittivityVacuum, lhs.getParameterisation());
            lhs.setDielectricPermittivityEM(dielectricPermittivityEMtemp);
        }
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
    auto size = conductivityEM.size();
    auto distEM = conductivityEM.getDistributionPtr();
    auto comm = distEM->getCommunicatorPtr();
    
    // create single distribution, only master process owns the complete vector (no distribution).
    
    int shotMaster = 0;
    auto singleDist = std::make_shared<dmemo::SingleDistribution>( size, comm, shotMaster );
    
    //redistribute vector to master process
    // (this may cause memory issues for big models)
    conductivityEM.redistribute(singleDist);    
    //reduce local array (size of local array is !=0 only for master process)
    commInterShot->sumArray(conductivityEM.getLocalValues());    
    //redistribute vector to former partition
    conductivityEM.redistribute(distEM);
        
    dielectricPermittivityEM.redistribute(singleDist);
    commInterShot->sumArray(dielectricPermittivityEM.getLocalValues());
    dielectricPermittivityEM.redistribute(distEM);
    
    porosity.redistribute(singleDist);
    commInterShot->sumArray(porosity.getLocalValues());
    porosity.redistribute(distEM);
    
    saturation.redistribute(singleDist);
    commInterShot->sumArray(saturation.getLocalValues());
    saturation.redistribute(distEM);
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
    auto distBig = dielectricPermittivityEM.getDistributionPtr();
    auto dist = gradientPerShot.getDielectricPermittivityEM().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> recoverMatrix = modelEM.getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd));
    recoverMatrix.assignTranspose(recoverMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = modelEM.getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd), boundaryWidth);
    scai::lama::SparseVector<ValueType> recoverVector;
    recoverVector = 1.0 - eraseVector;
    
    scai::lama::DenseVector<ValueType> temp;
    
    temp = recoverMatrix * gradientPerShot.getDielectricPermittivityEM(); //transform pershot into big model
    temp *= recoverVector;
    dielectricPermittivityEM *= eraseVector;
    dielectricPermittivityEM += temp; //take over the values
  
    temp = recoverMatrix * gradientPerShot.getConductivityEM(); //transform pershot into big model
    temp *= recoverVector;
    conductivityEM *= eraseVector;
    conductivityEM += temp; //take over the values
    
    temp = recoverMatrix * gradientPerShot.getPorosity(); //transform pershot into big model
    temp *= recoverVector;
    porosity *= eraseVector;
    porosity += temp; //take over the values
    
    temp = recoverMatrix * gradientPerShot.getSaturation(); //transform pershot into big model
    temp *= recoverVector;
    saturation *= eraseVector;
    saturation += temp; //take over the values
}

/*! \brief Function for scaling the gradients with the modelEM parameter 
 *
 \param modelEM Abstract modelEM.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::scale(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM, KITGPI::Configuration::Configuration configEM)
{    
    ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
    ValueType const ConductivityReference = modelEM.getConductivityReference();
    ValueType maxValue;      
    
    IndexType scaleGradient = configEM.get<IndexType>("scaleGradient");
    if (scaleGradient != 0) {
        if (workflowEM.getInvertForSigmaEM() && conductivityEM.maxNorm() != 0) {  
            if (scaleGradient == 1) {
                maxValue = modelEM.getConductivityEM().maxNorm();
            } else if (scaleGradient == 2) {
                maxValue = configEM.get<ValueType>("upperSigmaEMTh") - configEM.get<ValueType>("lowerSigmaEMTh");
            }
            this->applyParameterisation(maxValue, ConductivityReference, modelEM.getParameterisation());
            conductivityEM *= 1 / conductivityEM.maxNorm() * maxValue;
        }  
        
        if (workflowEM.getInvertForEpsilonEM() && dielectricPermittivityEM.maxNorm() != 0) {
            if (scaleGradient == 1) {
                maxValue = modelEM.getDielectricPermittivityEM().maxNorm();
            } else if (scaleGradient == 2) {
                maxValue = (configEM.get<ValueType>("upperEpsilonEMrTh") - configEM.get<ValueType>("lowerEpsilonEMrTh")) * DielectricPermittivityVacuum;
            }      
            this->applyParameterisation(maxValue, DielectricPermittivityVacuum, modelEM.getParameterisation()); 
            dielectricPermittivityEM *= 1 / dielectricPermittivityEM.maxNorm() * maxValue;
        }
        
        if (workflowEM.getInvertForPorosity() && porosity.maxNorm() != 0) {
            if (modelEM.getParameterisation() == 1 || modelEM.getParameterisation() == 2) {
                if (scaleGradient == 1) {
                    maxValue = modelEM.getPorosity().maxNorm();
                } else if (scaleGradient == 2) {
                    maxValue = configEM.get<ValueType>("upperPorosityTh") - configEM.get<ValueType>("lowerPorosityTh");
                }
            }        
            porosity *= 1 / porosity.maxNorm() * maxValue;
        }    
    
        if (workflowEM.getInvertForSaturation() && saturation.maxNorm() != 0) { 
            if (modelEM.getParameterisation() == 1 || modelEM.getParameterisation() == 2) {
                if (scaleGradient == 1) {
                    maxValue = modelEM.getSaturation().maxNorm();
                } else if (scaleGradient == 2) {
                    maxValue = configEM.get<ValueType>("upperSaturationTh") - configEM.get<ValueType>("lowerSaturationTh");
                }
            }              
            saturation *= 1 / saturation.maxNorm() * maxValue;        
        }
    }
}

/*! \brief Function for normalizing the gradientEM
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::normalize()
{
    if (this->getNormalizeGradient()) {
        ValueType gradientMax = conductivityEM.maxNorm();
        if (gradientMax != 0)
            conductivityEM *= 1 / gradientMax;
        gradientMax = dielectricPermittivityEM.maxNorm();
        if (gradientMax != 0)
            dielectricPermittivityEM *= 1 / gradientMax;
        gradientMax = porosity.maxNorm();
        if (gradientMax != 0)
            porosity *= 1 / gradientMax;
        gradientMax = saturation.maxNorm();
        if (gradientMax != 0)
            saturation *= 1 / gradientMax;
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
        conductivityEM = DT * correlatedWavefields.getXcorrSigmaEM(); 
        
        temp = modelEM.getConductivityEM();
        ValueType const ConductivityReference = modelEM.getConductivityReference();
        this->gradientParameterisation(conductivityEM, temp, ConductivityReference, modelEM.getParameterisation());
    } else {
        this->initParameterisation(conductivityEM, ctx, distEM, 0.0);
    }
    
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        dielectricPermittivityEM = gradEpsilonEM;
        
        temp = modelEM.getDielectricPermittivityEM();
        ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
        this->gradientParameterisation(dielectricPermittivityEM, temp, DielectricPermittivityVacuum, modelEM.getParameterisation());
    } else {
        this->initParameterisation(dielectricPermittivityEM, ctx, distEM, 0.0);
    }
    
    if (workflowEM.getInvertForPorosity()) {
        // porosity and water saturation        
        scai::lama::DenseVector<ValueType> conductivityDePorosity;
        scai::lama::DenseVector<ValueType> dielectricPermittiviyDePorosity;
     
        // Based on complex refractive index model (CRIM)    
        dielectricPermittiviyDePorosity = this->getDielectricPermittiviyDePorosity(modelEM);             
        
        porosity = dielectricPermittiviyDePorosity * dielectricPermittivityEM;
        if (modelEM.getParameterisation() == 2) {
            // Based on Archie equation
            conductivityDePorosity = this->getConductivityDePorosity(modelEM);    
            conductivityDePorosity *= conductivityEM;
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
        
        saturation = dielectricPermittiviyDeSaturation * dielectricPermittivityEM;
        if (modelEM.getParameterisation() == 2) {
            // Based on Archie equation
            conductivityDeSaturation = this->getConductivityDeSaturation(modelEM); 
            conductivityDeSaturation *= conductivityEM;
            saturation += conductivityDeSaturation; 
        }           
    } else {
        this->initParameterisation(saturation, ctx, distEM, 0.0);
    }  
}

/*! \brief calculate the stabilizing functional of each model parameter */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::calcStabilizingFunctionalGradient(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelPrioriEM, KITGPI::Configuration::Configuration configEM, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{    
    scai::lama::DenseVector<ValueType> conductivityEMPrioritemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEMPrioritemp;
    ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
    ValueType const ConductivityReference = modelEM.getConductivityReference();
    
    conductivityEM = modelEM.getConductivityEM();  
    dielectricPermittivityEM = modelEM.getDielectricPermittivityEM();     
    conductivityEMPrioritemp = modelPrioriEM.getConductivityEM();  
    dielectricPermittivityEMPrioritemp = modelPrioriEM.getDielectricPermittivityEM();

    this->applyParameterisation(conductivityEM, ConductivityReference, modelEM.getParameterisation());       
    this->applyParameterisation(dielectricPermittivityEM, DielectricPermittivityVacuum, modelEM.getParameterisation()); 
    this->applyParameterisation(conductivityEMPrioritemp, ConductivityReference, modelEM.getParameterisation());       
    this->applyParameterisation(dielectricPermittivityEMPrioritemp, DielectricPermittivityVacuum, modelEM.getParameterisation()); 
        
    scai::hmemo::ContextPtr ctx = dielectricPermittivityEM.getContextPtr();
    scai::dmemo::DistributionPtr distEM = dielectricPermittivityEM.getDistributionPtr();
    
    if (workflowEM.getInvertForEpsilonEM()) {
        dielectricPermittivityEMPrioritemp = dielectricPermittivityEM - dielectricPermittivityEMPrioritemp;
        dielectricPermittivityEM = this->calcStabilizingFunctionalGradientPerModel(dielectricPermittivityEMPrioritemp, configEM, dataMisfitEM);
    } else {
        this->initParameterisation(dielectricPermittivityEM, ctx, distEM, 0.0);
    }
    
    if (workflowEM.getInvertForSigmaEM()) { 
        conductivityEMPrioritemp = conductivityEM - conductivityEMPrioritemp;
        conductivityEM = this->calcStabilizingFunctionalGradientPerModel(conductivityEMPrioritemp, configEM, dataMisfitEM);
    } else {
        this->initParameterisation(conductivityEM, ctx, distEM, 0.0);
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
    
    sigmaEM_temp = this->getConductivityEM();
    epsilonEM_temp = this->getDielectricPermittivityEM();
    porosity_temp = this->getPorosity();
    saturation_temp = this->getSaturation();
    
    scai::IndexType NZ = configEM.get<IndexType>("NZ");
    if (NZ == 1) {
        scai::IndexType NY = configEM.get<IndexType>("NY");
        scai::IndexType NX = porosity_temp.size() / NY;
        scai::IndexType spatialFDorder = configEM.get<IndexType>("spatialFDorder");
        
        KITGPI::Common::applyMedianFilterTo2DVector(sigmaEM_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(epsilonEM_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(porosity_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(saturation_temp, NX, NY, spatialFDorder);
        
        this->setConductivityEM(sigmaEM_temp);    
        this->setDielectricPermittivityEM(epsilonEM_temp);
        this->setPorosity(porosity_temp);    
        this->setSaturation(saturation_temp);
    }
}

template class KITGPI::Gradient::EMEM<float>;
template class KITGPI::Gradient::EMEM<double>;
