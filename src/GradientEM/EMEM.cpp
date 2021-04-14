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
    conductivityEM *= rhs;
    dielectricPermittivityEM *= rhs;
    porosity *= rhs;
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
    conductivityEM *= rhs;
    dielectricPermittivityEM *= rhs;
    porosity *= rhs;
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
    scai::lama::DenseVector<ValueType> temp;  
    if (lhs.getParameterisation() == 2 || lhs.getParameterisation() == 1) { 
        temp = lhs.getPorosity() - rhs.getPorosity();
        lhs.setPorosity(temp);
        temp = lhs.getSaturation() - rhs.getSaturation();
        lhs.setSaturation(temp);         
    } else {
        scai::lama::DenseVector<ValueType> conductivityEMtemp;
        scai::lama::DenseVector<ValueType> dielectricPermittivityEMtemp;   
        ValueType const DielectricPermittivityVacuum = lhs.getDielectricPermittivityVacuum();
        ValueType const ConductivityReference = lhs.getConductivityReference();
        
        conductivityEMtemp = lhs.getConductivityEM();  
        dielectricPermittivityEMtemp = lhs.getDielectricPermittivityEM();
        
        this->exParameterisation(conductivityEMtemp, ConductivityReference, lhs.getParameterisation());       
        this->exParameterisation(dielectricPermittivityEMtemp, DielectricPermittivityVacuum, lhs.getParameterisation());  
                 
        conductivityEMtemp -= rhs.getConductivityEM();  
        dielectricPermittivityEMtemp -= rhs.getDielectricPermittivityEM(); 
        
        this->deParameterisation(conductivityEMtemp, ConductivityReference, lhs.getParameterisation());       
        this->deParameterisation(dielectricPermittivityEMtemp, DielectricPermittivityVacuum, lhs.getParameterisation());  
        
        lhs.setConductivityEM(conductivityEMtemp);
        lhs.setDielectricPermittivityEM(dielectricPermittivityEMtemp);
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
 \param gradientPerShot gradient per shot
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::sumGradientPerShot(KITGPI::Gradient::GradientEM<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth)
{
    auto distBig = dielectricPermittivityEM.getDistributionPtr();
    auto dist = gradientPerShot.getDielectricPermittivityEM().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix;
    scai::lama::DenseVector<ValueType> weightingVector(dist, 1.0);
    scai::lama::DenseVector<ValueType> weightingVectorBig(distBig, 0.0);
    IndexType numCuts = cutCoordinates.size();
    for (IndexType index=0; index < numCuts; index++) {        
        shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(index));
        shrinkMatrix.assignTranspose(shrinkMatrix);
        weightingVectorBig += shrinkMatrix * weightingVector;
    }
    weightingVectorBig = 1.0 / weightingVectorBig; // the weighting of overlapping area
    Common::replaceInvalid<ValueType>(weightingVectorBig, 0.0);
    
    shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd));
    shrinkMatrix.assignTranspose(shrinkMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = this->getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd), boundaryWidth);
    eraseVector *= weightingVectorBig;
    scai::lama::SparseVector<ValueType> restoreVector;
    restoreVector = 1.0 - eraseVector;
    
    scai::lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * gradientPerShot.getDielectricPermittivityEM(); //transform pershot into big model
    temp *= restoreVector;
    dielectricPermittivityEM *= eraseVector;
    dielectricPermittivityEM += temp; //take over the values
  
    temp = shrinkMatrix * gradientPerShot.getConductivityEM(); //transform pershot into big model
    temp *= restoreVector;
    conductivityEM *= eraseVector;
    conductivityEM += temp; //take over the values
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
    if (scaleGradient == 1) {
        if (workflowEM.getInvertForSigmaEM() && conductivityEM.maxNorm() != 0) {  
            scai::lama::DenseVector<ValueType> conductivityEMtemp;
            conductivityEMtemp = modelEM.getConductivityEM();  
            this->exParameterisation(conductivityEMtemp, ConductivityReference, modelEM.getParameterisation()); 
            conductivityEM *= 1 / conductivityEM.maxNorm() * conductivityEMtemp.maxNorm();
        }    
        if (workflowEM.getInvertForEpsilonEM() && dielectricPermittivityEM.maxNorm() != 0) {
            scai::lama::DenseVector<ValueType> dielectricPermittivityEMtemp;
            dielectricPermittivityEMtemp = modelEM.getDielectricPermittivityEM();        
            this->exParameterisation(dielectricPermittivityEMtemp, DielectricPermittivityVacuum, modelEM.getParameterisation()); 
            dielectricPermittivityEM *= 1 / dielectricPermittivityEM.maxNorm() * dielectricPermittivityEMtemp.maxNorm();
        }
        if (workflowEM.getInvertForPorosity() && porosity.maxNorm() != 0) {        
            porosity *= 1 / porosity.maxNorm() * modelEM.getPorosity().maxNorm();
        }
        if (workflowEM.getInvertForSaturation() && saturation.maxNorm() != 0) {       
            saturation *= 1 / saturation.maxNorm() * modelEM.getSaturation().maxNorm();        
        } 
    } else if (scaleGradient == 2) {
        if (workflowEM.getInvertForSigmaEM() && conductivityEM.maxNorm() != 0) { 
            maxValue = configEM.get<ValueType>("upperSigmaEMTh") - configEM.get<ValueType>("lowerSigmaEMTh");
            this->exParameterisation(maxValue, ConductivityReference, modelEM.getParameterisation());
            conductivityEM *= 1 / conductivityEM.maxNorm() * maxValue;
        }    
        if (workflowEM.getInvertForEpsilonEM() && dielectricPermittivityEM.maxNorm() != 0) {
            maxValue = configEM.get<ValueType>("upperEpsilonEMrTh") - configEM.get<ValueType>("lowerEpsilonEMrTh");
            this->exParameterisation(maxValue, DielectricPermittivityVacuum, modelEM.getParameterisation());
            dielectricPermittivityEM *= 1 / dielectricPermittivityEM.maxNorm() * maxValue;
        }
        if (workflowEM.getInvertForPorosity() && porosity.maxNorm() != 0) {  
            maxValue = configEM.get<ValueType>("upperPorosityTh") - configEM.get<ValueType>("lowerPorosityTh");
            porosity *= 1 / porosity.maxNorm() * maxValue; 
        }
        if (workflowEM.getInvertForSaturation() && saturation.maxNorm() != 0) {                 
            maxValue = configEM.get<ValueType>("upperSaturationTh") - configEM.get<ValueType>("lowerSaturationTh");
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

    this->exParameterisation(conductivityEM, ConductivityReference, modelEM.getParameterisation());       
    this->exParameterisation(dielectricPermittivityEM, DielectricPermittivityVacuum, modelEM.getParameterisation()); 
    this->exParameterisation(conductivityEMPrioritemp, ConductivityReference, modelEM.getParameterisation());       
    this->exParameterisation(dielectricPermittivityEMPrioritemp, DielectricPermittivityVacuum, modelEM.getParameterisation()); 
        
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
    
    bool useStreamConfig = configEM.get<bool>("useStreamConfig");
    scai::IndexType NX;
    scai::IndexType NY;
    if (!useStreamConfig) {
        NX = configEM.get<IndexType>("NX");
        NY = configEM.get<IndexType>("NY");
    } else {
        KITGPI::Configuration::Configuration configBigEM(configEM.get<std::string>("streamConfigFilename"));
        NX = configBigEM.get<IndexType>("NX");
        NY = configBigEM.get<IndexType>("NY");
    }          
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

template class KITGPI::Gradient::EMEM<float>;
template class KITGPI::Gradient::EMEM<double>;
