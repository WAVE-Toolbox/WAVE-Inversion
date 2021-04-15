#include "ViscoEMEM.hpp"
#include <scai/dmemo/SingleDistribution.hpp>

using namespace scai;

/*! \brief Constructor that is using the Configuration class
 *
 \param ctx Context for the Calculation
 \param distEM Distribution
 */
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType>::ViscoEMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM)
{
    equationTypeEM = "viscoemem";
    init(ctx, distEM, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

/*! \brief Initialisation with zeros
 *
 \param ctx Context for the Calculation
 \param distEM Distribution
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM)
{
    init(ctx, distEM, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

/*! \brief Initialisation that is generating a homogeneous modelEM
 *
 *  Generates a homogeneous modelEM, which will be initialized by the two given scalar values.
 \param ctx Context
 \param distEM Distribution
 \param conductivityEM_const conductivityEM given as Scalar
 \param dielectricPermittivityEM_const dielectricPermittivityEM given as Scalar
 \param tauConductivityEM_const tauConductivityEM given as Scalar
 \param tauDielectricPermittivityEM_const tauDielectricPermittivityEM given as Scalar
 \param porosity_const porosity given as Scalar
 \param saturation_const saturation given as Scalar
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, ValueType conductivityEM_const, ValueType dielectricPermittivityEM_const, ValueType tauConductivityEM_const, ValueType tauDielectricPermittivityEM_const, ValueType porosity_const, ValueType saturation_const)
{
    this->initParameterisation(conductivityEM, ctx, distEM, conductivityEM_const);
    this->initParameterisation(dielectricPermittivityEM, ctx, distEM, dielectricPermittivityEM_const);
    this->initParameterisation(tauConductivityEM, ctx, distEM, tauConductivityEM_const);
    this->initParameterisation(tauDielectricPermittivityEM, ctx, distEM, tauDielectricPermittivityEM_const);
    this->initParameterisation(porosity, ctx, distEM, porosity_const);
    this->initParameterisation(saturation, ctx, distEM, saturation_const);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::ViscoEMEM<ValueType>::ViscoEMEM(const ViscoEMEM &rhs)
{
    equationTypeEM = rhs.equationTypeEM;
    conductivityEM = rhs.conductivityEM;
    dielectricPermittivityEM = rhs.dielectricPermittivityEM;
    tauConductivityEM = rhs.tauConductivityEM;
    tauDielectricPermittivityEM = rhs.tauDielectricPermittivityEM;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
}

/*! \brief Write modelEM to an external file
 *
 \param filename For the tauConductivityEM ".tauSigmaEMr.mtx" and for tauDielectricPermittivityEM ".tauEpsilonEM.mtx" is added.
 \param fileFormat format of output file
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::write(std::string filename, IndexType fileFormat, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) const
{
    if (workflowEM.getInvertForSigmaEM()) {
        std::string filenameConductivityEM = filename + ".sigmaEM";
        this->writeParameterisation(conductivityEM, filenameConductivityEM, fileFormat);
    }
    
    if (workflowEM.getInvertForEpsilonEM()) {       
        std::string filenameEpsilonEM = filename + ".epsilonEMr";
        this->writeParameterisation(dielectricPermittivityEM, filenameEpsilonEM, fileFormat);
    }
    
    if (workflowEM.getInvertForTauSigmaEM()) {
        std::string filenameTauSigmaEM = filename + ".tauSigmaEMr";
        this->writeParameterisation(tauConductivityEM, filenameTauSigmaEM, fileFormat);
    }

    if (workflowEM.getInvertForTauEpsilonEM()) {
        std::string filenameTauEpsilonEM = filename + ".tauEpsilonEM";
        this->writeParameterisation(tauDielectricPermittivityEM, filenameTauEpsilonEM, fileFormat);
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

/*! \brief Get equationTypeEM (viscoemem)
 */
template <typename ValueType>
std::string KITGPI::Gradient::ViscoEMEM<ValueType>::getEquationType() const
{
    return (equationTypeEM);
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
    conductivityEM *= rhs;
    dielectricPermittivityEM *= rhs;
    tauConductivityEM *= rhs;
    tauDielectricPermittivityEM *= rhs;
    porosity *= rhs;
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
    conductivityEM += rhs.conductivityEM;
    dielectricPermittivityEM += rhs.dielectricPermittivityEM;
    tauConductivityEM += rhs.tauConductivityEM;
    tauDielectricPermittivityEM += rhs.tauDielectricPermittivityEM;
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
    conductivityEM -= rhs.conductivityEM;
    dielectricPermittivityEM -= rhs.dielectricPermittivityEM;
    tauConductivityEM -= rhs.tauConductivityEM;
    tauDielectricPermittivityEM -= rhs.tauDielectricPermittivityEM;
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
    conductivityEM = rhs.conductivityEM;
    dielectricPermittivityEM = rhs.dielectricPermittivityEM;
    tauConductivityEM = rhs.tauConductivityEM;
    tauDielectricPermittivityEM = rhs.tauDielectricPermittivityEM;
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
    conductivityEM = rhs.getConductivityEM();
    dielectricPermittivityEM = rhs.getDielectricPermittivityEM();
    tauConductivityEM = rhs.getTauConductivityEM();
    tauDielectricPermittivityEM = rhs.getTauDielectricPermittivityEM();
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

    conductivityEM -= rhs.getConductivityEM();
    dielectricPermittivityEM -= rhs.getDielectricPermittivityEM();
    tauConductivityEM -= rhs.getTauConductivityEM();
    tauDielectricPermittivityEM -= rhs.getTauDielectricPermittivityEM();
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

    conductivityEM += rhs.getConductivityEM();
    dielectricPermittivityEM += rhs.getDielectricPermittivityEM();
    tauConductivityEM += rhs.getTauConductivityEM();
    tauDielectricPermittivityEM += rhs.getTauDielectricPermittivityEM();
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
    conductivityEM *= rhs;
    dielectricPermittivityEM *= rhs;
    tauConductivityEM *= rhs;
    tauDielectricPermittivityEM *= rhs;
    porosity *= rhs;
    saturation *= rhs;
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradientEM which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::timesAssign(scai::lama::Vector<ValueType> const &rhs)
{
    conductivityEM *= rhs;
    dielectricPermittivityEM *= rhs;
    tauConductivityEM *= rhs;
    tauDielectricPermittivityEM *= rhs;
    porosity *= rhs;
    saturation *= rhs;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract modelEM.
 \param rhs Abstract gradientEM which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> &lhs, KITGPI::Gradient::GradientEM<ValueType> const &rhs)
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
        scai::lama::DenseVector<ValueType> tauConductivityEMtemp;
        scai::lama::DenseVector<ValueType> tauDielectricPermittivityEMtemp; 
        ValueType const DielectricPermittivityVacuum = lhs.getDielectricPermittivityVacuum();
        ValueType const ConductivityReference = lhs.getConductivityReference();    
        ValueType const TauDielectricPermittivityReference = lhs.getTauDielectricPermittivityReference();     
        ValueType const TauConductivityReference = lhs.getTauConductivityReference(); 
        
        conductivityEMtemp = lhs.getConductivityEM();  
        dielectricPermittivityEMtemp = lhs.getDielectricPermittivityEM(); 
        tauConductivityEMtemp = lhs.getTauConductivityEM();  
        tauDielectricPermittivityEMtemp = lhs.getTauDielectricPermittivityEM(); 

        this->exParameterisation(conductivityEMtemp, ConductivityReference, lhs.getParameterisation());       
        this->exParameterisation(dielectricPermittivityEMtemp, DielectricPermittivityVacuum, lhs.getParameterisation());  
        this->exParameterisation(tauConductivityEMtemp, TauConductivityReference, lhs.getParameterisation());       
        this->exParameterisation(tauDielectricPermittivityEMtemp, TauDielectricPermittivityReference, lhs.getParameterisation()); 
                 
        conductivityEMtemp -= rhs.getConductivityEM();  
        dielectricPermittivityEMtemp -= rhs.getDielectricPermittivityEM(); 
        tauConductivityEMtemp -= rhs.getTauConductivityEM();  
        tauDielectricPermittivityEMtemp -= rhs.getTauDielectricPermittivityEM(); 
        
        this->deParameterisation(conductivityEMtemp, ConductivityReference, lhs.getParameterisation());       
        this->deParameterisation(dielectricPermittivityEMtemp, DielectricPermittivityVacuum, lhs.getParameterisation());
        this->deParameterisation(tauConductivityEMtemp, TauConductivityReference, lhs.getParameterisation());       
        this->deParameterisation(tauDielectricPermittivityEMtemp, TauDielectricPermittivityReference, lhs.getParameterisation()); 
        
        lhs.setConductivityEM(conductivityEMtemp);
        lhs.setDielectricPermittivityEM(dielectricPermittivityEMtemp);
        lhs.setTauConductivityEM(tauConductivityEMtemp);
        lhs.setTauDielectricPermittivityEM(tauDielectricPermittivityEMtemp);
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
    auto size = dielectricPermittivityEM.size();
    auto distEM = dielectricPermittivityEM.getDistributionPtr();
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
    
    tauConductivityEM.redistribute(singleDist);
    commInterShot->sumArray(tauConductivityEM.getLocalValues());
    tauConductivityEM.redistribute(distEM);
    
    tauDielectricPermittivityEM.redistribute(singleDist);
    commInterShot->sumArray(tauDielectricPermittivityEM.getLocalValues());
    tauDielectricPermittivityEM.redistribute(distEM);
    
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
void KITGPI::Gradient::ViscoEMEM<ValueType>::sumGradientPerShot(KITGPI::Gradient::GradientEM<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth)
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
    
    temp = shrinkMatrix * gradientPerShot.getTauDielectricPermittivityEM(); //transform pershot into big model
    temp *= restoreVector;
    tauDielectricPermittivityEM *= eraseVector;
    tauDielectricPermittivityEM += temp; //take over the values
  
    temp = shrinkMatrix * gradientPerShot.getTauConductivityEM(); //transform pershot into big model
    temp *= restoreVector;
    tauConductivityEM *= eraseVector;
    tauConductivityEM += temp; //take over the values
}

/*! \brief Function for scaling the gradients with the modelEM parameter 
 *
 \param modelEM Abstract modelEM.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::scale(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM, KITGPI::Configuration::Configuration configEM)
{    
    ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum(); 
    ValueType const ConductivityReference = modelEM.getConductivityReference();    
    ValueType const TauDielectricPermittivityReference = modelEM.getTauDielectricPermittivityReference();      
    ValueType const TauConductivityReference = modelEM.getTauConductivityReference(); 
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
        if (workflowEM.getInvertForTauSigmaEM() && tauConductivityEM.maxNorm() != 0) {
            scai::lama::DenseVector<ValueType> tauConductivityEMtemp;
            tauConductivityEMtemp = modelEM.getTauConductivityEM();  
            this->exParameterisation(tauConductivityEMtemp, TauConductivityReference, modelEM.getParameterisation()); 
            tauConductivityEM *= 1 / tauConductivityEM.maxNorm() * tauConductivityEMtemp.maxNorm();
        }    
        if (workflowEM.getInvertForTauEpsilonEM() && tauDielectricPermittivityEM.maxNorm() != 0) {
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityEMtemp;
            tauDielectricPermittivityEMtemp = modelEM.getTauDielectricPermittivityEM(); 
            this->exParameterisation(tauDielectricPermittivityEMtemp, TauDielectricPermittivityReference, modelEM.getParameterisation()); 
            tauDielectricPermittivityEM *= 1 / tauDielectricPermittivityEM.maxNorm() * tauDielectricPermittivityEMtemp.maxNorm();
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
        if (workflowEM.getInvertForTauSigmaEM() && tauConductivityEM.maxNorm() != 0) {
            maxValue = configEM.get<ValueType>("upperTauSigmaEMrTh") - configEM.get<ValueType>("lowerTauSigmaEMrTh");  
            this->exParameterisation(maxValue, TauConductivityReference, modelEM.getParameterisation()); 
            tauConductivityEM *= 1 / tauConductivityEM.maxNorm() * maxValue;
        }    
        if (workflowEM.getInvertForTauEpsilonEM() && tauDielectricPermittivityEM.maxNorm() != 0) {
            maxValue = configEM.get<ValueType>("upperTauEpsilonEMTh") - configEM.get<ValueType>("lowerTauEpsilonEMTh");  
            this->exParameterisation(maxValue, TauDielectricPermittivityReference, modelEM.getParameterisation()); 
            tauDielectricPermittivityEM *= 1 / tauDielectricPermittivityEM.maxNorm() * maxValue;
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
void KITGPI::Gradient::ViscoEMEM<ValueType>::normalize()
{
    if (this->getNormalizeGradient()) {
        ValueType gradientMax = conductivityEM.maxNorm();
        if (gradientMax != 0)
            conductivityEM *= 1 / gradientMax;
        gradientMax = dielectricPermittivityEM.maxNorm();
        if (gradientMax != 0)
            dielectricPermittivityEM *= 1 / gradientMax;
        gradientMax = tauConductivityEM.maxNorm();
        if (gradientMax != 0)
            tauConductivityEM *= 1 / gradientMax;
        gradientMax = tauDielectricPermittivityEM.maxNorm();
        if (gradientMax != 0)
            tauDielectricPermittivityEM *= 1 / gradientMax;
        gradientMax = porosity.maxNorm();
        if (gradientMax != 0)
            porosity *= 1 / gradientMax;
        gradientMax = saturation.maxNorm();
        if (gradientMax != 0)
            saturation *= 1 / gradientMax;
    }
}

/*! \brief Function for calculating the viscoemem gradients from the cross correlation and the modelEM parameter 
 *
 \param modelEM Abstract modelEM.
 \param correlatedWavefields Abstract xCorr.
 \param DT Temporal discretization 
 \param workflowEM 
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
void KITGPI::Gradient::ViscoEMEM<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, ValueType DT, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{    
    scai::lama::DenseVector<ValueType> gradEpsilonEMoptical;
    scai::lama::DenseVector<ValueType> gradConductivityEMoptical;
    scai::lama::DenseVector<ValueType> gradEpsilonEMstatic;
    scai::lama::DenseVector<ValueType> gradTauEpsilonEM;
    scai::lama::DenseVector<ValueType> temp;     
    ValueType const TauDielectricPermittivityReference = modelEM.getTauDielectricPermittivityReference();   
    ValueType const TauConductivityReference = modelEM.getTauConductivityReference(); 
    
    gradEpsilonEMoptical = DT * correlatedWavefields.getXcorrEpsilonEM();
    
    scai::hmemo::ContextPtr ctx = gradEpsilonEMoptical.getContextPtr();
    scai::dmemo::DistributionPtr distEM = gradEpsilonEMoptical.getDistributionPtr();
//     std::cout << "estimateParameter modelEM.getNumRelaxationMechanisms() = " << modelEM.getNumRelaxationMechanisms() << "\n" << std::endl;
//     std::cout << "estimateParameter modelEM.getTauDisplacementEM() = " << modelEM.getTauDisplacementEM() << "\n" << std::endl;
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) { 
        gradConductivityEMoptical = DT * correlatedWavefields.getXcorrSigmaEM();  
    }
    
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {              
//         \nabla f_1(\varepsilon^s_e) =& -\sum_{l=1}^L \frac{L \tau_{Dl}}{{\varepsilon^s_e}^2 \tau_{\varepsilon e l}} \left( \tau_{Dl} r_{\varepsilon e l} + r_{\sigma e l} \right)
//         \nabla f_1(\tau_{\varepsilon e l}) =& - \frac{L \tau_{Dl}}{\varepsilon^s_e \tau_{\varepsilon e l}^2} \left( \tau_{Dl} r_{\varepsilon e l} + r_{\sigma e l} \right)
        temp = DT * correlatedWavefields.getXcorrREpsilonEM();    
        temp *= modelEM.getTauDisplacementEM();                      
        temp += DT * correlatedWavefields.getXcorrRSigmaEM();     
        temp *= modelEM.getTauDisplacementEM();            
        temp *= modelEM.getNumRelaxationMechanisms();      
        temp /= modelEM.getDielectricPermittivityEM();       
        temp /= modelEM.getTauDielectricPermittivityEM();  
        temp = - temp;    
        gradEpsilonEMstatic = temp / modelEM.getDielectricPermittivityEM();
        gradTauEpsilonEM = temp / modelEM.getTauDielectricPermittivityEM(); 
        // In case that tauDielectricPermittivityEM = 0 in air layer
        Common::replaceInvalid<ValueType>(gradEpsilonEMstatic, 0.0);
        Common::replaceInvalid<ValueType>(gradTauEpsilonEM, 0.0);
    }
    
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {      
        conductivityEM = gradConductivityEMoptical;
        temp = modelEM.getTauConductivityEM() * gradEpsilonEMoptical;
        conductivityEM += temp;  
        
        temp = modelEM.getConductivityEM();
        ValueType const ConductivityReference = modelEM.getConductivityReference();
        this->gradientParameterisation(conductivityEM, temp, ConductivityReference, modelEM.getParameterisation());
    } else {
        this->initParameterisation(conductivityEM, ctx, distEM, 0.0);
    }
    
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        dielectricPermittivityEM = gradEpsilonEMstatic; 
        temp = modelEM.getTauDielectricPermittivityEM() / modelEM.getTauDisplacementEM();               
        temp *= modelEM.getNumRelaxationMechanisms();      
        temp /= modelEM.getNumRelaxationMechanisms();
        temp *= gradConductivityEMoptical;  
        dielectricPermittivityEM += temp;
        temp = modelEM.getTauDielectricPermittivityEM() * modelEM.getNumRelaxationMechanisms();      
        temp /= modelEM.getNumRelaxationMechanisms();      
        temp = 1 - temp;     
        temp *= gradEpsilonEMoptical;
        dielectricPermittivityEM += temp;
        
        temp = modelEM.getDielectricPermittivityEM();
        ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
        this->gradientParameterisation(dielectricPermittivityEM, temp, DielectricPermittivityVacuum, modelEM.getParameterisation());
    } else {
        this->initParameterisation(dielectricPermittivityEM, ctx, distEM, 0.0);
    }

    if (workflowEM.getInvertForTauSigmaEM()) {
        tauConductivityEM = modelEM.getConductivityEM() * gradEpsilonEMoptical;
        
        temp = modelEM.getTauConductivityEM();
        this->gradientParameterisation(tauConductivityEM, temp, TauConductivityReference, modelEM.getParameterisation()); 
    } else {
        this->initParameterisation(tauConductivityEM, ctx, distEM, 0.0);
    }
    
    if (workflowEM.getInvertForTauEpsilonEM()) {
        tauDielectricPermittivityEM = gradTauEpsilonEM;
        temp = modelEM.getDielectricPermittivityEM() / modelEM.getNumRelaxationMechanisms();    
        temp /= modelEM.getTauDisplacementEM(); 
        temp *= gradConductivityEMoptical;  
        tauDielectricPermittivityEM += temp;  
        temp = modelEM.getDielectricPermittivityEM() / modelEM.getNumRelaxationMechanisms();
        temp = - temp;   
        temp *= gradEpsilonEMoptical; 
        tauDielectricPermittivityEM += temp;
        
        temp = modelEM.getTauDielectricPermittivityEM();
        this->gradientParameterisation(tauDielectricPermittivityEM, temp, TauDielectricPermittivityReference, modelEM.getParameterisation()); 
    } else {
        this->initParameterisation(tauDielectricPermittivityEM, ctx, distEM, 0.0);
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
void KITGPI::Gradient::ViscoEMEM<ValueType>::calcStabilizingFunctionalGradient(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelPrioriEM, KITGPI::Configuration::Configuration configEM, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
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
    
    this->initParameterisation(tauConductivityEM, ctx, distEM, 0.0);
    this->initParameterisation(tauDielectricPermittivityEM, ctx, distEM, 0.0);
}

/*! \brief Apply a median filter to filter the extrame value of the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::applyMedianFilter(KITGPI::Configuration::Configuration configEM)
{
    scai::lama::DenseVector<ValueType> sigmaEM_temp;
    scai::lama::DenseVector<ValueType> epsilonEM_temp;
    scai::lama::DenseVector<ValueType> tauSigmaEM_temp;
    scai::lama::DenseVector<ValueType> tauEpsilonEM_temp;
    scai::lama::DenseVector<ValueType> porosity_temp;
    scai::lama::DenseVector<ValueType> saturation_temp;
    
    sigmaEM_temp = this->getConductivityEM();
    epsilonEM_temp = this->getDielectricPermittivityEM();
    tauSigmaEM_temp = this->getTauConductivityEM();
    tauEpsilonEM_temp = this->getTauDielectricPermittivityEM();
    porosity_temp = this->getPorosity();
    saturation_temp = this->getSaturation();
    
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialFDorder = configEM.get<IndexType>("spatialFDorder");
    
    KITGPI::Common::applyMedianFilterTo2DVector(sigmaEM_temp, NX, NY, spatialFDorder);
    KITGPI::Common::applyMedianFilterTo2DVector(epsilonEM_temp, NX, NY, spatialFDorder);
    KITGPI::Common::applyMedianFilterTo2DVector(tauSigmaEM_temp, NX, NY, spatialFDorder);
    KITGPI::Common::applyMedianFilterTo2DVector(tauEpsilonEM_temp, NX, NY, spatialFDorder);
    KITGPI::Common::applyMedianFilterTo2DVector(porosity_temp, NX, NY, spatialFDorder);
    KITGPI::Common::applyMedianFilterTo2DVector(saturation_temp, NX, NY, spatialFDorder);
    
    this->setConductivityEM(sigmaEM_temp);    
    this->setDielectricPermittivityEM(epsilonEM_temp);
    this->setTauConductivityEM(tauSigmaEM_temp);    
    this->setTauDielectricPermittivityEM(tauEpsilonEM_temp);
    this->setPorosity(porosity_temp);    
    this->setSaturation(saturation_temp);
}

template class KITGPI::Gradient::ViscoEMEM<float>;
template class KITGPI::Gradient::ViscoEMEM<double>;
