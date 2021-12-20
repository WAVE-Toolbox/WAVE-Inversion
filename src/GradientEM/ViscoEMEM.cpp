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
    this->initParameterisation(reflectivity, ctx, dist, 0.0);
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
    reflectivity = rhs.reflectivity;
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
    this->resetParameter(reflectivity);
}

/*! \brief Write model to an external file
 *
 \param filename For the tauElectricConductivity ".tauSigmaEMr.mtx" and for tauDielectricPermittivity ".tauEpsilonEM.mtx" is added.
 \param fileFormat format of output file
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::write(std::string filename, IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const
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
    
    if (workflow.getInvertForReflectivity()) { 
        std::string filenameReflectivity = filename + ".reflectivity";
        this->writeParameterisation(reflectivity, filenameReflectivity, fileFormat);
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
    if (workflowInner.getInvertForReflectivity()) 
        reflectivity *= rhs;

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
    reflectivity += rhs.reflectivity;

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
    reflectivity -= rhs.reflectivity;

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
    reflectivity = rhs.reflectivity;

    return *this;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::assign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    electricConductivity = rhs.getElectricConductivity();
    dielectricPermittivity = rhs.getDielectricPermittivity();
    tauElectricConductivity = rhs.getTauElectricConductivity();
    tauDielectricPermittivity = rhs.getTauDielectricPermittivity();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    reflectivity = rhs.getReflectivity();
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    electricConductivity -= rhs.getElectricConductivity();
    dielectricPermittivity -= rhs.getDielectricPermittivity();
    tauElectricConductivity -= rhs.getTauElectricConductivity();
    tauDielectricPermittivity -= rhs.getTauDielectricPermittivity();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    reflectivity -= rhs.getReflectivity();
}

/*! \brief Function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{

    electricConductivity += rhs.getElectricConductivity();
    dielectricPermittivity += rhs.getDielectricPermittivity();
    tauElectricConductivity += rhs.getTauElectricConductivity();
    tauDielectricPermittivity += rhs.getTauDielectricPermittivity();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    reflectivity += rhs.getReflectivity();
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
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
    if (workflowInner.getInvertForReflectivity()) 
        reflectivity *= rhs;
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
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
    reflectivity *= rhs;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract model.
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs)
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
        scai::lama::DenseVector<ValueType> tauElectricConductivitytemp;
        scai::lama::DenseVector<ValueType> tauDielectricPermittivitytemp; 
        ValueType const DielectricPermittivityVacuum = lhs.getDielectricPermittivityVacuum();
        ValueType const ElectricConductivityReference = lhs.getElectricConductivityReference();    
        ValueType const TauDielectricPermittivityReference = lhs.getTauDielectricPermittivityReference();     
        ValueType const TauElectricConductivityReference = lhs.getTauElectricConductivityReference(); 
        
        if (workflowInner.getInvertForSigmaEM()) {       
            scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective; 
            if (lhs.getEffectiveParameterisation()) {
                dielectricPermittivityRealEffective = lhs.getDielectricPermittivityRealEffective();
                
                electricConductivitytemp = lhs.getElectricConductivityRealEffective();  
                this->applyParameterisation(electricConductivitytemp, ElectricConductivityReference, lhs.getParameterisation()); 
                electricConductivitytemp -= rhs.getElectricConductivity();    
                this->deleteParameterisation(electricConductivitytemp, ElectricConductivityReference, lhs.getParameterisation()); 
                
                IndexType calculateType = 2;
                scai::lama::DenseVector<ValueType> electricConductivityStatic = lhs.getElectricConductivityStatic(dielectricPermittivityRealEffective, electricConductivitytemp, calculateType); 
                lhs.setElectricConductivity(electricConductivityStatic);
            } else {
                if (!workflowInner.getInvertForEpsilonEM())
                    dielectricPermittivityRealEffective = lhs.getDielectricPermittivityRealEffective();
                
                electricConductivitytemp = lhs.getElectricConductivity();  
                this->applyParameterisation(electricConductivitytemp, ElectricConductivityReference, lhs.getParameterisation()); 
                electricConductivitytemp -= rhs.getElectricConductivity();    
                this->deleteParameterisation(electricConductivitytemp, ElectricConductivityReference, lhs.getParameterisation()); 
                lhs.setElectricConductivity(electricConductivitytemp);
                
                if (!workflowInner.getInvertForEpsilonEM()) {     
                    IndexType calculateType = 1;
                    scai::lama::DenseVector<ValueType> dielectricPermittivityStatic = lhs.getDielectricPermittivityStatic(dielectricPermittivityRealEffective, electricConductivitytemp, calculateType); 
                    lhs.setDielectricPermittivity(dielectricPermittivityStatic);
                }
            }
        }
        if (workflowInner.getInvertForEpsilonEM()) {
            scai::lama::DenseVector<ValueType> electricConductivityRealEffective; 
            if (lhs.getEffectiveParameterisation()) { 
                electricConductivityRealEffective = lhs.getElectricConductivityRealEffective();
                
                dielectricPermittivitytemp = lhs.getDielectricPermittivityRealEffective(); 
                this->applyParameterisation(dielectricPermittivitytemp, DielectricPermittivityVacuum, lhs.getParameterisation());  
                dielectricPermittivitytemp -= rhs.getDielectricPermittivity(); 
                this->deleteParameterisation(dielectricPermittivitytemp, DielectricPermittivityVacuum, lhs.getParameterisation());
                
                IndexType calculateType = 2;
                scai::lama::DenseVector<ValueType> dielectricPermittivityStatic = lhs.getDielectricPermittivityStatic(dielectricPermittivitytemp, electricConductivityRealEffective, calculateType); 
                lhs.setDielectricPermittivity(dielectricPermittivityStatic);
            } else {
                if (!workflowInner.getInvertForSigmaEM())
                    electricConductivityRealEffective = lhs.getElectricConductivityRealEffective();
                
                dielectricPermittivitytemp = lhs.getDielectricPermittivity(); 
                this->applyParameterisation(dielectricPermittivitytemp, DielectricPermittivityVacuum, lhs.getParameterisation());  
                dielectricPermittivitytemp -= rhs.getDielectricPermittivity(); 
                this->deleteParameterisation(dielectricPermittivitytemp, DielectricPermittivityVacuum, lhs.getParameterisation());
                lhs.setDielectricPermittivity(dielectricPermittivitytemp);
                
                if (!workflowInner.getInvertForSigmaEM()) {    
                    IndexType calculateType = 1;
                    scai::lama::DenseVector<ValueType> electricConductivityStatic = lhs.getElectricConductivityStatic(dielectricPermittivitytemp, electricConductivityRealEffective, calculateType); 
                    lhs.setElectricConductivity(electricConductivityStatic);
                }
            }
        }
        if (workflowInner.getInvertForTauSigmaEM()) {            
            tauElectricConductivitytemp = lhs.getTauElectricConductivity();  
            this->applyParameterisation(tauElectricConductivitytemp, TauElectricConductivityReference, lhs.getParameterisation());  
            tauElectricConductivitytemp -= rhs.getTauElectricConductivity();  
            this->deleteParameterisation(tauElectricConductivitytemp, TauElectricConductivityReference, lhs.getParameterisation()); 
            lhs.setTauElectricConductivity(tauElectricConductivitytemp);
        }
        if (workflowInner.getInvertForTauEpsilonEM()) {
            scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective; 
            scai::lama::DenseVector<ValueType> electricConductivityRealEffective;  
            if (!workflowInner.getInvertForEpsilonEM() && workflowInner.getInvertForSigmaEM()) {  
                dielectricPermittivityRealEffective = lhs.getDielectricPermittivityRealEffective();
            } else if (!workflowInner.getInvertForSigmaEM() && workflowInner.getInvertForEpsilonEM()) {
                electricConductivityRealEffective = lhs.getElectricConductivityRealEffective();
            }
            
            tauDielectricPermittivitytemp = lhs.getTauDielectricPermittivity();          
            this->applyParameterisation(tauDielectricPermittivitytemp, TauDielectricPermittivityReference, lhs.getParameterisation());                  
            tauDielectricPermittivitytemp -= rhs.getTauDielectricPermittivity();                    
            this->deleteParameterisation(tauDielectricPermittivitytemp, TauDielectricPermittivityReference, lhs.getParameterisation());         
            lhs.setTauDielectricPermittivity(tauDielectricPermittivitytemp);
                        
            if (!workflowInner.getInvertForEpsilonEM() && workflowInner.getInvertForSigmaEM()) {     
                IndexType calculateType = 2;
                electricConductivityRealEffective = lhs.getElectricConductivityRealEffective();
                scai::lama::DenseVector<ValueType> dielectricPermittivityStatic = lhs.getDielectricPermittivityStatic(dielectricPermittivityRealEffective, electricConductivityRealEffective, calculateType); 
                lhs.setDielectricPermittivity(dielectricPermittivityStatic);
            } else if (!workflowInner.getInvertForSigmaEM() && workflowInner.getInvertForEpsilonEM()) {
                IndexType calculateType = 2;
                dielectricPermittivityRealEffective = lhs.getDielectricPermittivityRealEffective();
                scai::lama::DenseVector<ValueType> electricConductivityStatic = lhs.getElectricConductivityStatic(dielectricPermittivityRealEffective, electricConductivityRealEffective, calculateType); 
                lhs.setElectricConductivity(electricConductivityStatic);
            }
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
    
    reflectivity.redistribute(singleDist);
    commInterShot->sumArray(reflectivity.getLocalValues());
    reflectivity.redistribute(dist);
}

/*! \brief If stream configuration is used, set a gradient per shot into the big gradient
 \param model model
 \param gradientPerShot gradient per shot
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::sumGradientPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd)
{
    auto distBig = dielectricPermittivity.getDistributionPtr();
    auto dist = gradientPerShot.getDielectricPermittivity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> recoverMatrix;
    
    recoverMatrix = model.getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd));
    recoverMatrix.assignTranspose(recoverMatrix);
        
    scai::lama::DenseVector<ValueType> temp;
    scai::lama::DenseVector<ValueType> weightingVector;
    scai::IndexType NY = modelCoordinates.getNY();
    
    if (workflowInner.getInvertForEpsilonEM()) {
        weightingVector = gradientPerShot.calcWeightingVector(gradientPerShot.getDielectricPermittivity(), NY, shotInd);
        temp = weightingVector * gradientPerShot.getDielectricPermittivity();  
        temp = recoverMatrix * temp;
        dielectricPermittivity += temp; 
    }
  
    if (workflowInner.getInvertForSigmaEM()) {
        weightingVector = gradientPerShot.calcWeightingVector(gradientPerShot.getElectricConductivity(), NY, shotInd);
        temp = weightingVector * gradientPerShot.getElectricConductivity();  
        temp = recoverMatrix * temp;
        electricConductivity += temp; 
    }
        
    if (workflowInner.getInvertForTauEpsilonEM()) {
        weightingVector = gradientPerShot.calcWeightingVector(gradientPerShot.getTauDielectricPermittivity(), NY, shotInd);
        temp = weightingVector * gradientPerShot.getTauDielectricPermittivity();  
        temp = recoverMatrix * temp;
        tauDielectricPermittivity += temp; 
    }
  
    if (workflowInner.getInvertForTauSigmaEM()) {
        weightingVector = gradientPerShot.calcWeightingVector(gradientPerShot.getTauElectricConductivity(), NY, shotInd);
        temp = weightingVector * gradientPerShot.getTauElectricConductivity();  
        temp = recoverMatrix * temp;
        tauElectricConductivity += temp; 
    }
    if (workflowInner.getInvertForPorosity()) {
        weightingVector = gradientPerShot.calcWeightingVector(gradientPerShot.getPorosity(), NY, shotInd);
        temp = weightingVector * gradientPerShot.getPorosity();  
        temp = recoverMatrix * temp;
        porosity += temp; 
    }
        
    if (workflowInner.getInvertForSaturation()) {
        weightingVector = gradientPerShot.calcWeightingVector(gradientPerShot.getSaturation(), NY, shotInd);
        temp = weightingVector * gradientPerShot.getSaturation();  
        temp = recoverMatrix * temp;
        saturation += temp; 
    }
        
    if (workflowInner.getInvertForReflectivity()) {
        weightingVector = gradientPerShot.calcWeightingVector(gradientPerShot.getReflectivity(), NY, shotInd);
        temp = weightingVector * gradientPerShot.getReflectivity();  
        temp = recoverMatrix * temp;
        reflectivity += temp; 
    }
}

/*! \brief Smooth gradient by Gaussian window
\param modelCoordinates coordinate class object of the model
*/
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::smooth(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config)
{
    if (config.get<IndexType>("NZ") == 1) {
        scai::IndexType smoothGradient = config.getAndCatch("smoothGradient", 0);
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity.size() / NY; // NX is different in stream configuration
        if (smoothGradient != 0) {
            HOST_PRINT(commAll, "\nApply Gaussian filter to gradient\n");
            scai::lama::DenseVector<ValueType> vector2Dpadded;
            if (workflowInner.getInvertForEpsilonEM()) {
                KITGPI::Common::pad2DVector(dielectricPermittivity, vector2Dpadded, NX, NY, PX, PY); 
                dielectricPermittivity = GaussianKernel * vector2Dpadded;
            }
            if (workflowInner.getInvertForSigmaEM()) {
                KITGPI::Common::pad2DVector(electricConductivity, vector2Dpadded, NX, NY, PX, PY); 
                electricConductivity = GaussianKernel * vector2Dpadded;
            }
            if (workflowInner.getInvertForTauEpsilonEM()) {
                KITGPI::Common::pad2DVector(tauDielectricPermittivity, vector2Dpadded, NX, NY, PX, PY); 
                tauDielectricPermittivity = GaussianKernel * vector2Dpadded;
            }
            if (workflowInner.getInvertForTauSigmaEM()) {
                KITGPI::Common::pad2DVector(tauElectricConductivity, vector2Dpadded, NX, NY, PX, PY); 
                tauElectricConductivity = GaussianKernel * vector2Dpadded;
            }
            if (workflowInner.getInvertForPorosity()) {
                KITGPI::Common::pad2DVector(porosity, vector2Dpadded, NX, NY, PX, PY); 
                porosity = GaussianKernel * vector2Dpadded;
            }
            if (workflowInner.getInvertForSaturation()) {
                KITGPI::Common::pad2DVector(saturation, vector2Dpadded, NX, NY, PX, PY); 
                saturation = GaussianKernel * vector2Dpadded;
            }
        }
    }
}

/*! \brief Smooth gradient by Gaussian window
\param modelCoordinates coordinate class object of the model
*/
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::applyMedianFilter(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config)
{
    if (config.get<IndexType>("NZ") == 1) {
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porosity.size() / NY; // NX is different in stream configuration
    
        HOST_PRINT(commAll, "\nApply median filter to gradient\n");
        scai::IndexType spatialLength = config.get<IndexType>("spatialFDorder");
        
        if (workflowInner.getInvertForEpsilonEM())
            KITGPI::Common::applyMedianFilterTo2DVector(dielectricPermittivity, NX, NY, spatialLength);
        if (workflowInner.getInvertForSigmaEM())
            KITGPI::Common::applyMedianFilterTo2DVector(electricConductivity, NX, NY, spatialLength);
        if (workflowInner.getInvertForTauEpsilonEM())
            KITGPI::Common::applyMedianFilterTo2DVector(tauDielectricPermittivity, NX, NY, spatialLength);
        if (workflowInner.getInvertForTauSigmaEM())
            KITGPI::Common::applyMedianFilterTo2DVector(tauElectricConductivity, NX, NY, spatialLength);
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
void KITGPI::Gradient::ViscoEMEM<ValueType>::scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{    
    ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum(); 
    ValueType const ElectricConductivityReference = model.getElectricConductivityReference();    
    ValueType const TauDielectricPermittivityReference = model.getTauDielectricPermittivityReference();      
    ValueType const TauElectricConductivityReference = model.getTauElectricConductivityReference(); 
    ValueType maxValue = 1;      
          
    IndexType scaleGradient = config.get<IndexType>("scaleGradient");
    if (scaleGradient != 0) {       
        if (workflow.getInvertForSigmaEM() && electricConductivity.maxNorm() != 0) {  
            if (scaleGradient == 1) {
                if (model.getEffectiveParameterisation()) {
                    maxValue = model.getElectricConductivityRealEffective().maxNorm();
                } else {
                    maxValue = model.getElectricConductivity().maxNorm();
                }
            } else if (scaleGradient == 2) {
                maxValue = config.get<ValueType>("upperSigmaEMTh") - config.get<ValueType>("lowerSigmaEMTh");
            }
            this->applyParameterisation(maxValue, ElectricConductivityReference, model.getParameterisation());
            electricConductivity *= 1 / electricConductivity.maxNorm() * maxValue;
        }  
        
        if (workflow.getInvertForEpsilonEM() && dielectricPermittivity.maxNorm() != 0) {
            if (scaleGradient == 1) {
                if (model.getEffectiveParameterisation()) {
                    maxValue = model.getDielectricPermittivityRealEffective().maxNorm();
                } else {
                    maxValue = model.getDielectricPermittivity().maxNorm();
                }
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
                ValueType relaxationTime_ref = 1.0 / (2.0 * M_PI * model.getCenterFrequencyCPML());
                maxValue = (config.get<ValueType>("upperTauSigmaEMrTh") - config.get<ValueType>("lowerTauSigmaEMrTh")) * relaxationTime_ref;
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
 */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::applyEnergyPreconditioning(ValueType epsilonHessian, scai::IndexType saveApproxHessian, std::string filename, scai::IndexType fileFormat)
{
    // see Nuber et al., 2015., Enhancement of near-surface elastic full waveform inversion results in regions of low sensitivities.
    scai::lama::DenseVector<ValueType> approxHessian;  
       
    if (workflowInner.getInvertForSigmaEM() && electricConductivity.maxNorm() != 0) {
        approxHessian = electricConductivity;
        approxHessian *= electricConductivity;
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, filename + ".sigmaEM", fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        electricConductivity *= approxHessian; 
    }
    
    if (workflowInner.getInvertForEpsilonEM() && dielectricPermittivity.maxNorm() != 0) {
        approxHessian = dielectricPermittivity;
        approxHessian *= dielectricPermittivity;
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, filename + ".epsilonEMr", fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        dielectricPermittivity *= approxHessian; 
    }
    
    if (workflowInner.getInvertForTauSigmaEM() && tauElectricConductivity.maxNorm() != 0) {
        approxHessian = tauElectricConductivity;
        approxHessian *= tauElectricConductivity;
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, filename + ".tauSigmaEMr", fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        tauElectricConductivity *= approxHessian; 
    }
    
    if (workflowInner.getInvertForTauEpsilonEM() && tauDielectricPermittivity.maxNorm() != 0) {
        approxHessian = tauDielectricPermittivity;
        approxHessian *= tauDielectricPermittivity;
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, filename + ".tauEpsilonEM", fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        tauDielectricPermittivity *= approxHessian; 
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
        gradientMax = reflectivity.maxNorm();
        if (gradientMax != 0)
            reflectivity *= 1 / gradientMax;
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
void KITGPI::Gradient::ViscoEMEM<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{    
    scai::lama::DenseVector<ValueType> gradEpsilonEffectiveOptical;
    scai::lama::DenseVector<ValueType> gradSigmaEffectiveOptical;
    scai::lama::DenseVector<ValueType> gradEpsilonStatic;
    scai::lama::DenseVector<ValueType> gradTauEpsilonEM;
    scai::lama::DenseVector<ValueType> gradEpsilonStatic1;
    scai::lama::DenseVector<ValueType> gradSigmaStatic;
    scai::lama::DenseVector<ValueType> gradTauEpsilonEM1;
    scai::lama::DenseVector<ValueType> gradTauSigmaEM;
    scai::lama::DenseVector<ValueType> gradEpsilonEffective;
    scai::lama::DenseVector<ValueType> gradSigmaEffective;
    scai::lama::DenseVector<ValueType> gradTauEpsilonEM2;
    scai::lama::DenseVector<ValueType> gradTauSigmaEM1;
    scai::lama::DenseVector<ValueType> temp; 
    scai::lama::DenseVector<ValueType> a_temp;
    scai::lama::DenseVector<ValueType> b_temp;  
    scai::lama::DenseVector<ValueType> c_temp;     
    ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum(); 
    ValueType const ElectricConductivityReference = model.getElectricConductivityReference();    
    ValueType const TauDielectricPermittivityReference = model.getTauDielectricPermittivityReference();      
    ValueType const TauElectricConductivityReference = model.getTauElectricConductivityReference();  
    ValueType w_average = 0;
    IndexType numRelaxationMechanisms = model.getNumRelaxationMechanisms();
    std::vector<ValueType> relaxationTime;          // = 1 / ( 2 * Pi * f_relax )
    ValueType a_average = 0;
    ValueType b_average = 0;
    ValueType w_ref = 2.0 * M_PI * model.getCenterFrequencyCPML();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        relaxationTime.push_back(1.0 / (2.0 * M_PI * model.getRelaxationFrequency()[l])); // = 1 / ( 2 * Pi * f_relax )
        w_average += 1.0 / relaxationTime[l];
        a_average += (w_ref * w_ref * relaxationTime[l] * relaxationTime[l] / (1 + w_ref * w_ref * relaxationTime[l] * relaxationTime[l]));
        b_average += (w_ref * w_ref * relaxationTime[l] / (1 + w_ref * w_ref * relaxationTime[l] * relaxationTime[l]));
    }
    w_average /= numRelaxationMechanisms;
    a_average /= numRelaxationMechanisms;
    b_average /= numRelaxationMechanisms;
    a_temp = a_average * model.getTauDielectricPermittivity();
    a_temp = 1 - a_temp;
    b_temp = b_average * model.getTauDielectricPermittivity();
    c_temp = b_average * model.getTauElectricConductivity();
    c_temp = a_temp - c_temp;
    c_temp = 1.0 / c_temp;

    gradEpsilonEffectiveOptical = DT * correlatedWavefields.getXcorrEpsilonEM();
    
    scai::hmemo::ContextPtr ctx = gradEpsilonEffectiveOptical.getContextPtr();
    scai::dmemo::DistributionPtr dist = gradEpsilonEffectiveOptical.getDistributionPtr();
    
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForTauSigmaEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) { 
        // basic gradients
        gradSigmaEffectiveOptical = DT * correlatedWavefields.getXcorrSigmaEM(); 
                
        scai::lama::DenseVector<ValueType> gradREpsilonSigmaEM;   
        gradREpsilonSigmaEM = DT * correlatedWavefields.getXcorrREpsilonSigmaEM(); 
        gradREpsilonSigmaEM *= -numRelaxationMechanisms;      
        gradREpsilonSigmaEM /= model.getDielectricPermittivity();       
        gradREpsilonSigmaEM /= model.getTauDielectricPermittivity();  
        gradEpsilonStatic = gradREpsilonSigmaEM / model.getDielectricPermittivity();
        gradTauEpsilonEM = gradREpsilonSigmaEM / model.getTauDielectricPermittivity(); 
        
        // In case that tauDielectricPermittivity = 0 in air layer
        Common::replaceInvalid<ValueType>(gradEpsilonStatic, 0.0);
        Common::replaceInvalid<ValueType>(gradTauEpsilonEM, 0.0);
             
        // convert the gradient from (gradEpsilonEffectiveOptical, gradSigmaEffectiveOptical, gradEpsilonStatic, gradTauEpsilonEM) to (gradEpsilonStatic, gradSigmaStatic, gradTauEpsilonEM1, gradTauSigmaEM).
        gradSigmaStatic = gradSigmaEffectiveOptical;
        temp = model.getTauElectricConductivity() * gradEpsilonEffectiveOptical;
        gradSigmaStatic += temp;  
        
        gradEpsilonStatic1 = gradEpsilonStatic;         
        temp = gradSigmaEffectiveOptical * w_average;
        temp *= model.getTauDielectricPermittivity();  
        gradEpsilonStatic1 += temp;
        temp = 1 - model.getTauDielectricPermittivity();     
        temp *= gradEpsilonEffectiveOptical;
        gradEpsilonStatic1 += temp;
        
        gradTauSigmaEM = model.getElectricConductivity() * gradEpsilonEffectiveOptical;
        
        gradTauEpsilonEM1 = gradTauEpsilonEM;        
        temp = gradSigmaEffectiveOptical * w_average;
        temp *= model.getDielectricPermittivity();
        gradTauEpsilonEM1 += temp;  
        temp = gradEpsilonEffectiveOptical * model.getDielectricPermittivity(); 
        gradTauEpsilonEM1 -= temp;        
    }
    
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) { 
        if (model.getEffectiveParameterisation()) {
            scai::lama::DenseVector<ValueType> epsilonStatic1DeSigmaEffective;
            scai::lama::DenseVector<ValueType> sigmaStaticDeSigmaEffective;
            epsilonStatic1DeSigmaEffective = -c_temp * model.getTauElectricConductivity();
            sigmaStaticDeSigmaEffective = b_temp * epsilonStatic1DeSigmaEffective;
            sigmaStaticDeSigmaEffective = 1 - sigmaStaticDeSigmaEffective;
            
            electricConductivity = epsilonStatic1DeSigmaEffective * gradEpsilonStatic1;
            temp = sigmaStaticDeSigmaEffective * gradSigmaStatic;
            electricConductivity += temp;
        } else {
            electricConductivity = gradSigmaStatic;  
        }
        
        this->gradientParameterisation(electricConductivity, model.getElectricConductivity(), ElectricConductivityReference, model.getParameterisation());
    } else {
        this->initParameterisation(electricConductivity, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {        
        if (model.getEffectiveParameterisation()) {
            scai::lama::DenseVector<ValueType> epsilonStatic1DeEpsilonEffective;
            scai::lama::DenseVector<ValueType> sigmaStaticDeEpsilonEffective;
            epsilonStatic1DeEpsilonEffective = c_temp;
            sigmaStaticDeEpsilonEffective = -b_temp * epsilonStatic1DeEpsilonEffective;
            
            dielectricPermittivity = epsilonStatic1DeEpsilonEffective * gradEpsilonStatic1;
            temp = sigmaStaticDeEpsilonEffective * gradSigmaStatic;
            dielectricPermittivity += temp;
        } else {
            dielectricPermittivity = gradEpsilonStatic1;
        }
        
        this->gradientParameterisation(dielectricPermittivity, model.getDielectricPermittivity(), DielectricPermittivityVacuum, model.getParameterisation());
    } else {
        this->initParameterisation(dielectricPermittivity, ctx, dist, 0.0);
    }

    if (workflow.getInvertForTauSigmaEM()) {              
        if (model.getEffectiveParameterisation()) {
            scai::lama::DenseVector<ValueType> epsilonStatic1DeTauSigma1;
            scai::lama::DenseVector<ValueType> sigmaStaticDeTauSigma1;
            epsilonStatic1DeTauSigma1 = b_temp * model.getDielectricPermittivityRealEffective();
            temp = a_temp * model.getElectricConductivityRealEffective();
            epsilonStatic1DeTauSigma1 -= temp;
            epsilonStatic1DeTauSigma1 *= c_temp;
            epsilonStatic1DeTauSigma1 *= c_temp;
            sigmaStaticDeTauSigma1 = -b_temp * epsilonStatic1DeTauSigma1;
            
            tauElectricConductivity = epsilonStatic1DeTauSigma1 * gradEpsilonStatic1;
            temp = sigmaStaticDeTauSigma1 * gradSigmaStatic;
            tauElectricConductivity += temp;
            tauElectricConductivity += gradTauSigmaEM;
        } else {
            tauElectricConductivity = gradTauSigmaEM;
        }
        
        this->gradientParameterisation(tauElectricConductivity, model.getTauElectricConductivity(), TauElectricConductivityReference, model.getParameterisation()); 
    } else {
        this->initParameterisation(tauElectricConductivity, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForTauEpsilonEM()) {                 
        if (model.getEffectiveParameterisation()) {
            scai::lama::DenseVector<ValueType> epsilonStatic1DeTauEpsilon2;
            scai::lama::DenseVector<ValueType> sigmaStaticDeTauEpsilon2;
            epsilonStatic1DeTauEpsilon2 = model.getElectricConductivityRealEffective() * model.getTauElectricConductivity();
            epsilonStatic1DeTauEpsilon2 = model.getDielectricPermittivityRealEffective()- epsilonStatic1DeTauEpsilon2;
            temp = b_average * model.getTauElectricConductivity();
            temp += a_average;
            epsilonStatic1DeTauEpsilon2 *= temp;
            epsilonStatic1DeTauEpsilon2 *= c_temp;
            epsilonStatic1DeTauEpsilon2 *= c_temp;
            sigmaStaticDeTauEpsilon2 = -b_temp * epsilonStatic1DeTauEpsilon2;
            
            tauDielectricPermittivity = epsilonStatic1DeTauEpsilon2 * gradEpsilonStatic1;
            temp = sigmaStaticDeTauEpsilon2 * gradSigmaStatic;
            tauDielectricPermittivity += temp;
            tauDielectricPermittivity += gradTauEpsilonEM1;
        } else {
            tauDielectricPermittivity = gradTauEpsilonEM1;
        }
        
        this->gradientParameterisation(tauDielectricPermittivity, model.getTauDielectricPermittivity(), TauDielectricPermittivityReference, model.getParameterisation()); 
    } else {
        this->initParameterisation(tauDielectricPermittivity, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForPorosity()) {
        // porosity and water saturation        
        scai::lama::DenseVector<ValueType> conductivityDePorosity;
        scai::lama::DenseVector<ValueType> dielectricPermittiviyDePorosity;
     
        // Based on complex refractive index model (CRIM)    
        dielectricPermittiviyDePorosity = this->getDielectricPermittiviyDePorosity(model);             
        
        porosity = dielectricPermittiviyDePorosity * gradEpsilonEffectiveOptical;
        if (model.getParameterisation() == 2) {
            // Based on Archie equation
            conductivityDePorosity = this->getElectricConductivityDePorosity(model);    
            conductivityDePorosity *= gradSigmaEffectiveOptical;
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
        
        saturation = dielectricPermittiviyDeSaturation * gradEpsilonEffectiveOptical;
        if (model.getParameterisation() == 2) {
            // Based on Archie equation
            conductivityDeSaturation = this->getElectricConductivityDeSaturation(model); 
            conductivityDeSaturation *= gradSigmaEffectiveOptical;
            saturation += conductivityDeSaturation; 
        }               
    } else {
        this->initParameterisation(saturation, ctx, dist, 0.0);
    }  
    
    if (workflow.getInvertForReflectivity()) {
        reflectivity = -gradEpsilonEffectiveOptical;
    } else {
        this->initParameterisation(reflectivity, ctx, dist, 0.0);
    } 
}

/*! \brief calculate the stabilizing functional of each model parameter */
template <typename ValueType>
void KITGPI::Gradient::ViscoEMEM<ValueType>::calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPrioriEM, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Workflow::Workflow<ValueType> const &workflow)
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

template class KITGPI::Gradient::ViscoEMEM<float>;
template class KITGPI::Gradient::ViscoEMEM<double>;
