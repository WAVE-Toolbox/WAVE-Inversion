#include "ViscoEMEM.hpp"
#include <scai/dmemo/SingleDistribution.hpp>

using namespace scai;

/*! \brief Constructor that is using the Configuration class
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType>::ViscoEMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "viscoemem";
    init(ctx, dist, 0.0, 0.0, 0.0, 0.0);
}

/*! \brief Initialisation with zeros
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist, 0.0, 0.0, 0.0, 0.0);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param electricConductivity_const electricConductivity given as Scalar
 \param dielectricPermittivity_const dielectricPermittivity given as Scalar
 \param tauElectricConductivity_const tauElectricConductivity given as Scalar
 \param tauDielectricPermittivity_const tauDielectricPermittivity given as Scalar
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType electricConductivity_const, ValueType dielectricPermittivity_const, ValueType tauElectricConductivity_const, ValueType tauDielectricPermittivity_const)
{
    this->initParameterisation(electricConductivity, ctx, dist, electricConductivity_const);
    this->initParameterisation(dielectricPermittivity, ctx, dist, dielectricPermittivity_const);
    this->initParameterisation(tauElectricConductivity, ctx, dist, tauElectricConductivity_const);
    this->initParameterisation(tauDielectricPermittivity, ctx, dist, tauDielectricPermittivity_const);
    this->initParameterisation(porosity, ctx, dist, 0.0);
    this->initParameterisation(saturation, ctx, dist, 0.0);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType>::ViscoEMEM(const ViscoEMEM &rhs)
{
    equationType = rhs.equationType;
    electricConductivity = rhs.electricConductivity;
    dielectricPermittivity = rhs.dielectricPermittivity;
    tauElectricConductivity = rhs.tauElectricConductivity;
    tauDielectricPermittivity = rhs.tauDielectricPermittivity;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
}

/*! \brief Set all parameter to zero.
*/
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::resetGradient()
{
    this->resetParameter(electricConductivity);
    this->resetParameter(dielectricPermittivity);
    this->resetParameter(tauElectricConductivity);
    this->resetParameter(tauDielectricPermittivity);
    this->resetParameter(porosity);
    this->resetParameter(saturation);
}

/*! \brief Write model to an external file
 *
 \param filename For the tauElectricConductivity ".tauSigmaEMr.mtx" and for tauDielectricPermittivity ".tauEpsilonEM.mtx" is added.
 \param fileFormat format of output file
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::write(std::string filename, IndexType fileFormat, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) const
{
    if (workflow.getInvertForSigmaEM()) {
        std::string filenameElectricConductivity = filename + ".sigmaEM";
        this->writeParameterisation(electricConductivity, filenameElectricConductivity, fileFormat);
    }
    
    if (workflow.getInvertForEpsilonEM()) {       
        std::string filenameDielectricPermittivity = filename + ".epsilonEMr";
        this->writeParameterisation(dielectricPermittivity, filenameDielectricPermittivity, fileFormat);
    }
    
    if (workflow.getInvertForTauSigmaEM()) {
        std::string filenameTauSigmaEM = filename + ".tauSigmaEMr";
        this->writeParameterisation(tauElectricConductivity, filenameTauSigmaEM, fileFormat);
    }

    if (workflow.getInvertForTauEpsilonEM()) {
        std::string filenameTauEpsilonEM = filename + ".tauEpsilonEM";
        this->writeParameterisation(tauDielectricPermittivity, filenameTauEpsilonEM, fileFormat);
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

/*! \brief Get equationType (viscoemem)
 */
template <typename ValueType>
std::string KITGPI::Gradient::ViscoEMEM<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType> KITGPI::Gradient::ViscoEMEM<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Gradient::ViscoEMEM<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief Free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType> operator*(ValueType lhs, KITGPI::Gradient::ViscoEMEM<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType> &KITGPI::Gradient::ViscoEMEM<ValueType>::operator*=(ValueType const &rhs)
{
    if (workflowInner.getInvertForSigmaEM()) 
        electricConductivity *= rhs;
    if (workflowInner.getInvertForEpsilonEM()) 
        dielectricPermittivity *= rhs;
    if (workflowInner.getInvertForTauSigmaEM()) 
        tauElectricConductivity *= rhs;
    if (workflowInner.getInvertForTauEpsilonEM()) 
        tauDielectricPermittivity *= rhs;
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
KITGPI::Gradient::ViscoEMEM<ValueType> KITGPI::Gradient::ViscoEMEM<ValueType>::operator+(KITGPI::Gradient::ViscoEMEM<ValueType> const &rhs)
{
    KITGPI::Gradient::ViscoEMEM<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType> &KITGPI::Gradient::ViscoEMEM<ValueType>::operator+=(KITGPI::Gradient::ViscoEMEM<ValueType> const &rhs)
{
    electricConductivity += rhs.electricConductivity;
    dielectricPermittivity += rhs.dielectricPermittivity;
    tauElectricConductivity += rhs.tauElectricConductivity;
    tauDielectricPermittivity += rhs.tauDielectricPermittivity;
    porosity += rhs.porosity;
    saturation += rhs.saturation;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType> KITGPI::Gradient::ViscoEMEM<ValueType>::operator-(KITGPI::Gradient::ViscoEMEM<ValueType> const &rhs)
{
    KITGPI::Gradient::ViscoEMEM<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType> &KITGPI::Gradient::ViscoEMEM<ValueType>::operator-=(KITGPI::Gradient::ViscoEMEM<ValueType> const &rhs)
{
    electricConductivity -= rhs.electricConductivity;
    dielectricPermittivity -= rhs.dielectricPermittivity;
    tauElectricConductivity -= rhs.tauElectricConductivity;
    tauDielectricPermittivity -= rhs.tauDielectricPermittivity;
    porosity -= rhs.porosity;
    saturation -= rhs.saturation;

    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Gradient which is copied.
 */
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType> &KITGPI::Gradient::ViscoEMEM<ValueType>::operator=(KITGPI::Gradient::ViscoEMEM<ValueType> const &rhs)
{
    electricConductivity = rhs.electricConductivity;
    dielectricPermittivity = rhs.dielectricPermittivity;
    tauElectricConductivity = rhs.tauElectricConductivity;
    tauDielectricPermittivity = rhs.tauDielectricPermittivity;
    porosity = rhs.porosity;
    saturation = rhs.saturation;

    return *this;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::assign(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    electricConductivity = rhs.getElectricConductivity();
    dielectricPermittivity = rhs.getDielectricPermittivity();
    tauElectricConductivity = rhs.getTauElectricConductivity();
    tauDielectricPermittivity = rhs.getTauDielectricPermittivity();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::minusAssign(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    electricConductivity -= rhs.getElectricConductivity();
    dielectricPermittivity -= rhs.getDielectricPermittivity();
    tauElectricConductivity -= rhs.getTauElectricConductivity();
    tauDielectricPermittivity -= rhs.getTauDielectricPermittivity();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
}

/*! \brief Function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::plusAssign(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{

    electricConductivity += rhs.getElectricConductivity();
    dielectricPermittivity += rhs.getDielectricPermittivity();
    tauElectricConductivity += rhs.getTauElectricConductivity();
    tauDielectricPermittivity += rhs.getTauDielectricPermittivity();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::timesAssign(ValueType const &rhs)
{
    if (workflowInner.getInvertForSigmaEM()) 
        electricConductivity *= rhs;
    if (workflowInner.getInvertForEpsilonEM()) 
        dielectricPermittivity *= rhs;
    if (workflowInner.getInvertForTauSigmaEM()) 
        tauElectricConductivity *= rhs;
    if (workflowInner.getInvertForTauEpsilonEM()) 
        tauDielectricPermittivity *= rhs;
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
void KITGPI::Gradient::ViscoEMEM<ValueType>::timesAssign(scai::lama::Vector<ValueType> const &rhs)
{
    electricConductivity *= rhs;
    dielectricPermittivity *= rhs;
    tauElectricConductivity *= rhs;
    tauDielectricPermittivity *= rhs;
    porosity *= rhs;
    saturation *= rhs;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract model.
 \param rhs Abstract gradientEM which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> &lhs, KITGPI::Gradient::GradientEM<ValueType> const &rhs)
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
        scai::lama::DenseVector<ValueType> electricConductivitytemp;
        scai::lama::DenseVector<ValueType> dielectricPermittivitytemp; 
        scai::lama::DenseVector<ValueType> tauElectricConductivitytemp;
        scai::lama::DenseVector<ValueType> tauDielectricPermittivitytemp; 
        ValueType const DielectricPermittivityVacuum = lhs.getDielectricPermittivityVacuum();
        ValueType const ElectricConductivityReference = lhs.getElectricConductivityReference();    
        ValueType const TauDielectricPermittivityReference = lhs.getTauDielectricPermittivityReference();     
        ValueType const TauElectricConductivityReference = lhs.getTauElectricConductivityReference(); 
        
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
        if (workflowInner.getInvertForTauSigmaEM()) {
            tauElectricConductivitytemp = lhs.getTauElectricConductivity();  
            this->applyParameterisation(tauElectricConductivitytemp, TauElectricConductivityReference, lhs.getParameterisation());  
            tauElectricConductivitytemp -= rhs.getTauElectricConductivity();  
            this->deleteParameterisation(tauElectricConductivitytemp, TauElectricConductivityReference, lhs.getParameterisation()); 
            lhs.setTauElectricConductivity(tauElectricConductivitytemp);
        }
        if (workflowInner.getInvertForTauEpsilonEM()) {
            tauDielectricPermittivitytemp = lhs.getTauDielectricPermittivity();          
            this->applyParameterisation(tauDielectricPermittivitytemp, TauDielectricPermittivityReference, lhs.getParameterisation());                  
            tauDielectricPermittivitytemp -= rhs.getTauDielectricPermittivity();                    
            this->deleteParameterisation(tauDielectricPermittivitytemp, TauDielectricPermittivityReference, lhs.getParameterisation());         
            lhs.setTauDielectricPermittivity(tauDielectricPermittivitytemp);
        }
    }
};

/*! \brief Function for summing the gradients of all shot domains
 *
 \param commInterShot inter shot communication pointer
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot)
{
    /*reduction between shot domains.
    each shot domain may have a different distribution of (gradient) vectors. This happens if geographer is used (different result for dist on each shot domain even for homogenous architecture) or on heterogenous architecture. In this case even the number of processes on each domain can vary. Therfore it is necessary that only one process per shot domain communicates all data.
    */
    
    //get information from distributed vector
    auto size = dielectricPermittivity.size();
    auto dist = dielectricPermittivity.getDistributionPtr();
    auto comm = dist->getCommunicatorPtr();
    
    // create single distribution, only master process owns the complete vector (no distribution).
    
    int shotMaster = 0;
    auto singleDist = std::make_shared<dmemo::SingleDistribution>( size, comm, shotMaster );
    
    //redistribute vector to master process
    // (this may cause memory issues for big models)
    electricConductivity.redistribute(singleDist);    
    //reduce local array (size of local array is !=0 only for master process)
    commInterShot->sumArray(electricConductivity.getLocalValues());    
    //redistribute vector to former partition
    electricConductivity.redistribute(dist);
        
    dielectricPermittivity.redistribute(singleDist);
    commInterShot->sumArray(dielectricPermittivity.getLocalValues());
    dielectricPermittivity.redistribute(dist);
    
    tauElectricConductivity.redistribute(singleDist);
    commInterShot->sumArray(tauElectricConductivity.getLocalValues());
    tauElectricConductivity.redistribute(dist);
    
    tauDielectricPermittivity.redistribute(singleDist);
    commInterShot->sumArray(tauDielectricPermittivity.getLocalValues());
    tauDielectricPermittivity.redistribute(dist);
    
    porosity.redistribute(singleDist);
    commInterShot->sumArray(porosity.getLocalValues());
    porosity.redistribute(dist);
    
    saturation.redistribute(singleDist);
    commInterShot->sumArray(saturation.getLocalValues());
    saturation.redistribute(dist);
}

/*! \brief If stream configuration is used, set a gradient per shot into the big gradient
 \param model model
 \param gradientPerShot gradient per shot
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::sumGradientPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &model, KITGPI::Gradient::GradientEM<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth)
{
    auto distBig = dielectricPermittivity.getDistributionPtr();
    auto dist = gradientPerShot.getDielectricPermittivity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> recoverMatrix;
    
    recoverMatrix = model.getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd));
    recoverMatrix.assignTranspose(recoverMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = model.getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd), boundaryWidth);
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
    
    temp = recoverMatrix * gradientPerShot.getTauDielectricPermittivity(); //transform pershot into big model
    temp *= recoverVector;
    tauDielectricPermittivity *= eraseVector;
    tauDielectricPermittivity += temp; //take over the values
  
    temp = recoverMatrix * gradientPerShot.getTauElectricConductivity(); //transform pershot into big model
    temp *= recoverVector;
    tauElectricConductivity *= eraseVector;
    tauElectricConductivity += temp; //take over the values
    
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
void KITGPI::Gradient::ViscoEMEM<ValueType>::scale(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{    
    ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum(); 
    ValueType const ElectricConductivityReference = model.getElectricConductivityReference();    
    ValueType const TauDielectricPermittivityReference = model.getTauDielectricPermittivityReference();      
    ValueType const TauElectricConductivityReference = model.getTauElectricConductivityReference(); 
    ValueType maxValue = 1;      
          
    IndexType scaleGradient = config.get<IndexType>("scaleGradient");
    if (workflow.getInvertForSigmaEM() && electricConductivity.maxNorm() != 0) {  
        if (scaleGradient == 1) {
            maxValue = model.getElectricConductivity().maxNorm();
        } else if (scaleGradient == 2) {
            maxValue = config.get<ValueType>("upperSigmaEMTh") - config.get<ValueType>("lowerSigmaEMTh");
        }
        this->applyParameterisation(maxValue, ElectricConductivityReference, model.getParameterisation());
        electricConductivity *= 1 / electricConductivity.maxNorm() * maxValue;
    }  
    
    if (workflow.getInvertForEpsilonEM() && dielectricPermittivity.maxNorm() != 0) {
        if (scaleGradient == 1) {
            maxValue = model.getDielectricPermittivity().maxNorm();
        } else if (scaleGradient == 2) {
            maxValue = (config.get<ValueType>("upperEpsilonEMrTh") - config.get<ValueType>("lowerEpsilonEMrTh")) * DielectricPermittivityVacuum;
        }      
        this->applyParameterisation(maxValue, DielectricPermittivityVacuum, model.getParameterisation()); 
        dielectricPermittivity *= 1 / dielectricPermittivity.maxNorm() * maxValue;
    }
    
    if (workflow.getInvertForTauSigmaEM() && tauElectricConductivity.maxNorm() != 0) {
        if (scaleGradient == 1) {
            maxValue = model.getTauElectricConductivity().maxNorm();
        } else if (scaleGradient == 2) {
            maxValue = (config.get<ValueType>("upperTauSigmaEMrTh") - config.get<ValueType>("lowerTauSigmaEMrTh")) * model.getTauElectricDisplacement();
        }
        this->applyParameterisation(maxValue, TauElectricConductivityReference, model.getParameterisation());
        tauElectricConductivity *= 1 / tauElectricConductivity.maxNorm() * maxValue;
    }    
    
    if (workflow.getInvertForTauEpsilonEM() && tauDielectricPermittivity.maxNorm() != 0) {
        if (scaleGradient == 1) {
            maxValue = model.getTauDielectricPermittivity().maxNorm();
        } else if (scaleGradient == 2) {
            maxValue = config.get<ValueType>("upperTauEpsilonEMTh") - config.get<ValueType>("lowerTauEpsilonEMTh");
        }      
        this->applyParameterisation(maxValue, TauDielectricPermittivityReference, model.getParameterisation()); 
        tauDielectricPermittivity *= 1 / tauDielectricPermittivity.maxNorm() * maxValue;
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

/*! \brief Function for normalizing the gradientEM
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::normalize()
{
    if (this->getNormalizeGradient()) {
        ValueType gradientMax = electricConductivity.maxNorm();
        if (gradientMax != 0)
            electricConductivity *= 1 / gradientMax;
        gradientMax = dielectricPermittivity.maxNorm();
        if (gradientMax != 0)
            dielectricPermittivity *= 1 / gradientMax;
        gradientMax = tauElectricConductivity.maxNorm();
        if (gradientMax != 0)
            tauElectricConductivity *= 1 / gradientMax;
        gradientMax = tauDielectricPermittivity.maxNorm();
        if (gradientMax != 0)
            tauDielectricPermittivity *= 1 / gradientMax;
        gradientMax = porosity.maxNorm();
        if (gradientMax != 0)
            porosity *= 1 / gradientMax;
        gradientMax = saturation.maxNorm();
        if (gradientMax != 0)
            saturation *= 1 / gradientMax;
    }
}

/*! \brief Function for calculating the viscoemem gradients from the cross correlation and the model parameter 
 *
 \param model Abstract model.
 \param correlatedWavefields Abstract xCorr.
 \param DT Temporal discretization 
 \param workflow 
 *
    \begin{equation}
    \label{eqn:GradientAdjoint12}
    \begin{align*}
    \nabla f_1(\varepsilon^{s'}_e) =& \pdv{\varepsilon^\infty_e}{\varepsilon^s_e} \nabla f_1(\varepsilon^\infty_e) + \pdv{\sigma^\infty_e}{\varepsilon^s_e} \nabla f_1(\sigma^\infty_e) + \nabla f_1(\varepsilon^s_e)\\
    \nabla f_1(\sigma^s_e) =& \pdv{\varepsilon^\infty_e}{\sigma^s_e} \nabla f_1(\varepsilon^\infty_e) + \nabla f_1(\sigma^\infty_e) \\
    \nabla f_1(\tau^{'}_{\varepsilon e l}) =& \pdv{\varepsilon^\infty_e}{\tau_{\varepsilon e l}} \nabla f_1(\varepsilon^\infty_e) + \pdv{\sigma^\infty_e}{\tau_{\varepsilon e l}} \nabla f_1(\sigma^\infty_e) + \nabla f_1(\tau_{\varepsilon e l})\\
    \nabla f_1(\tau_{\sigma e}) =& \pdv{\varepsilon^\infty_e}{\tau_{\sigma e}} \nabla f_1(\varepsilon^\infty_e) 
    \end{align*}
    \end{equation}
    where 
    \begin{equation}
    \begin{align*}
    \pdv{\varepsilon^\infty_e}{\varepsilon^s_e} =& \left( 1-\frac{1 }{L} \sum_{l=1}^L \tau_{\varepsilon e l} \right),\quad
    \pdv{\sigma^\infty_e}{\varepsilon^s_e} = \frac{1 }{L} \sum_{l=1}^L \frac{\tau_{\varepsilon e l}}{\tau_{Dl}}\\
    \pdv{\varepsilon^\infty_e}{\sigma^s_e} =& \tau_{\sigma e}\\
    \pdv{\varepsilon^\infty_e}{\tau_{\varepsilon e l}} =&  -\frac{\varepsilon^s_{e}}{L},\quad
    \pdv{\sigma^\infty_e}{\tau_{\varepsilon e l}} = \frac{\varepsilon^s_{e}}{L \tau_{Dl}}\\
    \pdv{\varepsilon^\infty_e}{\tau_{\sigma e}} =& \sigma^s_e
    \end{align*}
    \end{equation}
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model, ValueType DT, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow)
{    
    scai::lama::DenseVector<ValueType> gradEpsilonEMoptical;
    scai::lama::DenseVector<ValueType> gradElectricConductivityoptical;
    scai::lama::DenseVector<ValueType> gradEpsilonEMstatic;
    scai::lama::DenseVector<ValueType> gradTauEpsilonEM;
    scai::lama::DenseVector<ValueType> temp;     
    ValueType const TauDielectricPermittivityReference = model.getTauDielectricPermittivityReference();   
    ValueType const TauElectricConductivityReference = model.getTauElectricConductivityReference(); 
    
    gradEpsilonEMoptical = DT * correlatedWavefields.getXcorrEpsilonEM();
    
    scai::hmemo::ContextPtr ctx = gradEpsilonEMoptical.getContextPtr();
    scai::dmemo::DistributionPtr dist = gradEpsilonEMoptical.getDistributionPtr();
//     std::cout << "estimateParameter model.getNumRelaxationMechanisms() = " << model.getNumRelaxationMechanisms() << "\n" << std::endl;
//     std::cout << "estimateParameter model.getTauElectricDisplacement() = " << model.getTauElectricDisplacement() << "\n" << std::endl;
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) { 
        gradElectricConductivityoptical = DT * correlatedWavefields.getXcorrSigmaEM();  
    }
    
    if (workflow.getInvertForEpsilonEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {              
//         \nabla f_1(\varepsilon^s_e) =& -\sum_{l=1}^L \frac{L \tau_{Dl}}{{\varepsilon^s_e}^2 \tau_{\varepsilon e l}} \left( \tau_{Dl} r_{\varepsilon e l} + r_{\sigma e l} \right)
//         \nabla f_1(\tau_{\varepsilon e l}) =& - \frac{L \tau_{Dl}}{\varepsilon^s_e \tau_{\varepsilon e l}^2} \left( \tau_{Dl} r_{\varepsilon e l} + r_{\sigma e l} \right)
        temp = DT * correlatedWavefields.getXcorrREpsilonEM();    
        temp *= model.getTauElectricDisplacement();                      
        temp += DT * correlatedWavefields.getXcorrRSigmaEM();     
        temp *= model.getTauElectricDisplacement();            
        temp *= model.getNumRelaxationMechanisms();      
        temp /= model.getDielectricPermittivity();       
        temp /= model.getTauDielectricPermittivity();  
        temp = - temp;    
        gradEpsilonEMstatic = temp / model.getDielectricPermittivity();
        gradTauEpsilonEM = temp / model.getTauDielectricPermittivity(); 
        // In case that tauDielectricPermittivity = 0 in air layer
        Common::replaceInvalid<ValueType>(gradEpsilonEMstatic, 0.0);
        Common::replaceInvalid<ValueType>(gradTauEpsilonEM, 0.0);
    }
    
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {      
        electricConductivity = gradElectricConductivityoptical;
        temp = model.getTauElectricConductivity() * gradEpsilonEMoptical;
        electricConductivity += temp;  
        
        temp = model.getElectricConductivity();
        ValueType const ElectricConductivityReference = model.getElectricConductivityReference();
        this->gradientParameterisation(electricConductivity, temp, ElectricConductivityReference, model.getParameterisation());
    } else {
        this->initParameterisation(electricConductivity, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        dielectricPermittivity = gradEpsilonEMstatic; 
        temp = model.getTauDielectricPermittivity() / model.getTauElectricDisplacement();               
        temp *= model.getNumRelaxationMechanisms();      
        temp /= model.getNumRelaxationMechanisms();
        temp *= gradElectricConductivityoptical;  
        dielectricPermittivity += temp;
        temp = model.getTauDielectricPermittivity() * model.getNumRelaxationMechanisms();      
        temp /= model.getNumRelaxationMechanisms();      
        temp = 1 - temp;     
        temp *= gradEpsilonEMoptical;
        dielectricPermittivity += temp;
        
        temp = model.getDielectricPermittivity();
        ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum();
        this->gradientParameterisation(dielectricPermittivity, temp, DielectricPermittivityVacuum, model.getParameterisation());
    } else {
        this->initParameterisation(dielectricPermittivity, ctx, dist, 0.0);
    }

    if (workflow.getInvertForTauSigmaEM()) {
        tauElectricConductivity = model.getElectricConductivity() * gradEpsilonEMoptical;
        
        temp = model.getTauElectricConductivity();
        this->gradientParameterisation(tauElectricConductivity, temp, TauElectricConductivityReference, model.getParameterisation()); 
    } else {
        this->initParameterisation(tauElectricConductivity, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForTauEpsilonEM()) {
        tauDielectricPermittivity = gradTauEpsilonEM;
        temp = model.getDielectricPermittivity() / model.getNumRelaxationMechanisms();    
        temp /= model.getTauElectricDisplacement(); 
        temp *= gradElectricConductivityoptical;  
        tauDielectricPermittivity += temp;  
        temp = model.getDielectricPermittivity();
        temp = - temp;   
        temp *= gradEpsilonEMoptical; 
        tauDielectricPermittivity += temp;
        
        temp = model.getTauDielectricPermittivity();
        this->gradientParameterisation(tauDielectricPermittivity, temp, TauDielectricPermittivityReference, model.getParameterisation()); 
    } else {
        this->initParameterisation(tauDielectricPermittivity, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForPorosity()) {
        // porosity and water saturation        
        scai::lama::DenseVector<ValueType> conductivityDePorosity;
        scai::lama::DenseVector<ValueType> dielectricPermittiviyDePorosity;
     
        // Based on complex refractive index model (CRIM)    
        dielectricPermittiviyDePorosity = this->getDielectricPermittiviyDePorosity(model);             
        
        porosity = dielectricPermittiviyDePorosity * dielectricPermittivity;
        if (model.getParameterisation() == 2) {
            // Based on Archie equation
            conductivityDePorosity = this->getElectricConductivityDePorosity(model);    
            conductivityDePorosity *= electricConductivity;
            porosity += conductivityDePorosity; 
        }        
    } else {
        this->initParameterisation(porosity, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForSaturation()) {
        // porosity and water saturation        
        scai::lama::DenseVector<ValueType> conductivityDeSaturation;
        scai::lama::DenseVector<ValueType> dielectricPermittiviyDeSaturation;
               
        // Based on complex refractive index model (CRIM)            
        dielectricPermittiviyDeSaturation = this->getDielectricPermittiviyDeSaturation(model);             
        
        saturation = dielectricPermittiviyDeSaturation * dielectricPermittivity;
        if (model.getParameterisation() == 2) {
            // Based on Archie equation
            conductivityDeSaturation = this->getElectricConductivityDeSaturation(model); 
            conductivityDeSaturation *= electricConductivity;
            saturation += conductivityDeSaturation; 
        }               
    } else {
        this->initParameterisation(saturation, ctx, dist, 0.0);
    }  
}

/*! \brief calculate the stabilizing functional of each model parameter */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::calcStabilizingFunctionalGradient(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelPrioriEM, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow)
{    
    scai::lama::DenseVector<ValueType> electricConductivityPrioritemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivityPrioritemp;
    ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum();
    ValueType const ElectricConductivityReference = model.getElectricConductivityReference();
    
    electricConductivity = model.getElectricConductivity();  
    dielectricPermittivity = model.getDielectricPermittivity();     
    electricConductivityPrioritemp = modelPrioriEM.getElectricConductivity();  
    dielectricPermittivityPrioritemp = modelPrioriEM.getDielectricPermittivity();

    this->applyParameterisation(electricConductivity, ElectricConductivityReference, model.getParameterisation());       
    this->applyParameterisation(dielectricPermittivity, DielectricPermittivityVacuum, model.getParameterisation()); 
    this->applyParameterisation(electricConductivityPrioritemp, ElectricConductivityReference, model.getParameterisation());       
    this->applyParameterisation(dielectricPermittivityPrioritemp, DielectricPermittivityVacuum, model.getParameterisation()); 
        
    scai::hmemo::ContextPtr ctx = dielectricPermittivity.getContextPtr();
    scai::dmemo::DistributionPtr dist = dielectricPermittivity.getDistributionPtr();
    
    if (workflow.getInvertForEpsilonEM()) {
        dielectricPermittivityPrioritemp = dielectricPermittivity - dielectricPermittivityPrioritemp;
        dielectricPermittivity = this->calcStabilizingFunctionalGradientPerModel(dielectricPermittivityPrioritemp, config, dataMisfitEM);
    } else {
        this->initParameterisation(dielectricPermittivity, ctx, dist, 0.0);
    }
        
    if (workflow.getInvertForSigmaEM()) { 
        electricConductivityPrioritemp = electricConductivity - electricConductivityPrioritemp;
        electricConductivity = this->calcStabilizingFunctionalGradientPerModel(electricConductivityPrioritemp, config, dataMisfitEM);
    } else {
        this->initParameterisation(electricConductivity, ctx, dist, 0.0);
    }    
    
    this->initParameterisation(tauElectricConductivity, ctx, dist, 0.0);
    this->initParameterisation(tauDielectricPermittivity, ctx, dist, 0.0);
}

/*! \brief Apply a median filter to filter the extrame value of the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::applyMedianFilter(KITGPI::Configuration::Configuration config)
{
    scai::lama::DenseVector<ValueType> sigmaEM_temp;
    scai::lama::DenseVector<ValueType> epsilonEM_temp;
    scai::lama::DenseVector<ValueType> tauSigmaEM_temp;
    scai::lama::DenseVector<ValueType> tauEpsilonEM_temp;
    scai::lama::DenseVector<ValueType> porosity_temp;
    scai::lama::DenseVector<ValueType> saturation_temp;
    
    sigmaEM_temp = this->getElectricConductivity();
    epsilonEM_temp = this->getDielectricPermittivity();
    tauSigmaEM_temp = this->getTauElectricConductivity();
    tauEpsilonEM_temp = this->getTauDielectricPermittivity();
    porosity_temp = this->getPorosity();
    saturation_temp = this->getSaturation();
    
    scai::IndexType NZ = config.get<IndexType>("NZ");
    if (NZ == 1) {
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity_temp.size() / NY;
        scai::IndexType spatialFDorder = config.get<IndexType>("spatialFDorder");
            
        KITGPI::Common::applyMedianFilterTo2DVector(sigmaEM_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(epsilonEM_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(tauSigmaEM_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(tauEpsilonEM_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(porosity_temp, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(saturation_temp, NX, NY, spatialFDorder);
        
        this->setElectricConductivity(sigmaEM_temp);    
        this->setDielectricPermittivity(epsilonEM_temp);
        this->setTauElectricConductivity(tauSigmaEM_temp);    
        this->setTauDielectricPermittivity(tauEpsilonEM_temp);
        this->setPorosity(porosity_temp);    
        this->setSaturation(saturation_temp);
    }
}

template class KITGPI::Gradient::ViscoEMEM<float>;
template class KITGPI::Gradient::ViscoEMEM<double>;
