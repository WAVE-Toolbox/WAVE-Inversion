#include "EMEM.hpp"
#include <scai/dmemo/SingleDistribution.hpp>

using namespace scai;

/*! \brief Constructor that is using the Configuration class
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType>::EMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "emem";
    init(ctx, dist, 0.0, 0.0);
}

/*! \brief Initialisation with zeros
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist, 0.0, 0.0);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param electricConductivity_const electricConductivity given as Scalar
 \param dielectricPermittivity_const dielectricPermittivity given as Scalar
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType electricConductivity_const, ValueType dielectricPermittivity_const)
{
    this->initParameterisation(electricConductivity, ctx, dist, electricConductivity_const);
    this->initParameterisation(dielectricPermittivity, ctx, dist, dielectricPermittivity_const);
    this->initParameterisation(porosity, ctx, dist, 0.0);
    this->initParameterisation(saturation, ctx, dist, 0.0);
    this->initParameterisation(reflectivity, ctx, dist, 0.0);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::EMEM<ValueType>::EMEM(const EMEM &rhs)
{
    equationType = rhs.equationType;
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

/*! \brief Write model to an external file
 *
 \param filename For the electricConductivity ".sigmaEM.mtx" and for dielectricPermittivity ".epsilonEMr.mtx" is added.
 \param fileFormat format of output file
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::write(std::string filename, IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const
{
    if (workflow.getInvertForSigmaEM()) {
        std::string filenameElectricConductivity = filename + ".sigmaEM";
        this->writeParameterisation(electricConductivity, filenameElectricConductivity, fileFormat);
    }
    
    if (workflow.getInvertForEpsilonEM()) {
        std::string filenameDielectricPermittivity = filename + ".epsilonEMr";
        this->writeParameterisation(dielectricPermittivity, filenameDielectricPermittivity, fileFormat);
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

/*! \brief Get equationType (emem)
 */
template <typename ValueType>
std::string KITGPI::Gradient::EMEM<ValueType>::getEquationType() const
{
    return (equationType);
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
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::assign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    electricConductivity = rhs.getElectricConductivity();
    dielectricPermittivity = rhs.getDielectricPermittivity();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    reflectivity = rhs.getReflectivity();
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    electricConductivity -= rhs.getElectricConductivity();
    dielectricPermittivity -= rhs.getDielectricPermittivity();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    reflectivity -= rhs.getReflectivity();
}

/*! \brief Function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    electricConductivity += rhs.getElectricConductivity();
    dielectricPermittivity += rhs.getDielectricPermittivity();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    reflectivity += rhs.getReflectivity();
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
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
 \param rhs Abstract gradient which is subtracted.
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
 \param lhs Abstract model.
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs)
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
    auto dist = electricConductivity.getDistributionPtr();
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
void KITGPI::Gradient::EMEM<ValueType>::sumGradientPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth)
{
    auto distBig = dielectricPermittivity.getDistributionPtr();
    auto dist = gradientPerShot.getDielectricPermittivity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> recoverMatrix = model.getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd));
    recoverMatrix.assignTranspose(recoverMatrix);
    
//     scai::lama::SparseVector<ValueType> eraseVector = model.getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotInd), boundaryWidth);
//     scai::lama::SparseVector<ValueType> recoverVector;
//     recoverVector = 1.0 - eraseVector;
    
    scai::lama::DenseVector<ValueType> temp;
    
    temp = recoverMatrix * gradientPerShot.getDielectricPermittivity(); //transform pershot into big model
//     temp *= recoverVector;
//     dielectricPermittivity *= eraseVector;
    dielectricPermittivity += temp; //take over the values
  
    temp = recoverMatrix * gradientPerShot.getElectricConductivity(); //transform pershot into big model
//     temp *= recoverVector;
//     electricConductivity *= eraseVector;
    electricConductivity += temp; //take over the values
    
    temp = recoverMatrix * gradientPerShot.getPorosity(); //transform pershot into big model
//     temp *= recoverVector;
//     porosity *= eraseVector;
    porosity += temp; //take over the values
    
    temp = recoverMatrix * gradientPerShot.getSaturation(); //transform pershot into big model
//     temp *= recoverVector;
//     saturation *= eraseVector;
    saturation += temp; //take over the values
    
    temp = recoverMatrix * gradientPerShot.getReflectivity(); //transform pershot into big model
//     temp *= recoverVector;
//     reflectivity *= eraseVector;
    reflectivity += temp; //take over the values
}

/*! \brief Smooth gradient by Gaussian window
\param modelCoordinates coordinate class object of the model
*/
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::smooth(scai::dmemo::CommunicatorPtr commAll, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType FCmax)
{
    if (smoothGradient == 2 || smoothGradient == 3) {
        HOST_PRINT(commAll, "\nApply Gaussian smooth to gradient\n");
        scai::lama::DenseVector<ValueType> vector2Dpadded;
        if (workflowInner.getInvertForEpsilonEM()) {
            KITGPI::Common::pad2DVector(dielectricPermittivity, vector2Dpadded, ksize, modelCoordinates); 
            dielectricPermittivity = GaussianKernel * vector2Dpadded;
        }
        if (workflowInner.getInvertForSigmaEM()) {
            KITGPI::Common::pad2DVector(electricConductivity, vector2Dpadded, ksize, modelCoordinates); 
            electricConductivity = GaussianKernel * vector2Dpadded;
        }
        if (workflowInner.getInvertForPorosity()) {
            KITGPI::Common::pad2DVector(porosity, vector2Dpadded, ksize, modelCoordinates); 
            porosity = GaussianKernel * vector2Dpadded;
        }
        if (workflowInner.getInvertForSaturation()) {
            KITGPI::Common::pad2DVector(saturation, vector2Dpadded, ksize, modelCoordinates); 
            saturation = GaussianKernel * vector2Dpadded;
        }
    }
}

/*! \brief Function for scaling the gradients with the model parameter 
 *
 \param model Abstract model.
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{    
    ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum();
    ValueType const ElectricConductivityReference = model.getElectricConductivityReference();
    ValueType maxValue = 1;      
    
    IndexType scaleGradient = config.get<IndexType>("scaleGradient");
    if (scaleGradient != 0) {
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
void KITGPI::Gradient::EMEM<ValueType>::applyEnergyPreconditioning(ValueType epsilonHessian, scai::IndexType saveApproxHessian, std::string filename, scai::IndexType fileFormat)
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

/*! \brief Function for calculating the emem gradients from the cross correlation and the model parameter 
 *
 \param model Abstract model.
 \param correlatedWavefields Abstract xCorr.
 \param DT Temporal discretization 
 \param workflow 
 *
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{    
    scai::lama::DenseVector<ValueType> gradEpsilonEM;
    scai::lama::DenseVector<ValueType> gradSigmaEM;
    scai::lama::DenseVector<ValueType> temp;     
    
    gradEpsilonEM = DT * correlatedWavefields.getXcorrEpsilonEM();
    
    scai::hmemo::ContextPtr ctx = gradEpsilonEM.getContextPtr();
    scai::dmemo::DistributionPtr dist = gradEpsilonEM.getDistributionPtr();
    
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) { 
        gradSigmaEM = DT * correlatedWavefields.getXcorrSigmaEM(); 
        electricConductivity = gradSigmaEM;
        
        this->gradientParameterisation(electricConductivity, model.getElectricConductivity(), model.getElectricConductivityReference(), model.getParameterisation());
    } else {
        this->initParameterisation(electricConductivity, ctx, dist, 0.0);
    }
    
    if (workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation() || workflow.getInvertForReflectivity()) {
        dielectricPermittivity = gradEpsilonEM;
        
        this->gradientParameterisation(dielectricPermittivity, model.getDielectricPermittivity(), model.getDielectricPermittivityVacuum(), model.getParameterisation());
    } else {
        this->initParameterisation(dielectricPermittivity, ctx, dist, 0.0);
    }
        
    if (workflow.getInvertForPorosity()) {
        // porosity and water saturation        
        scai::lama::DenseVector<ValueType> conductivityDePorosity;
        scai::lama::DenseVector<ValueType> dielectricPermittiviyDePorosity;
     
        // Based on complex refractive index model (CRIM)    
        dielectricPermittiviyDePorosity = this->getDielectricPermittiviyDePorosity(model);             
        
        porosity = dielectricPermittiviyDePorosity * gradEpsilonEM;
        if (model.getParameterisation() == 2) {
            // Based on Archie equation
            conductivityDePorosity = this->getElectricConductivityDePorosity(model);    
            conductivityDePorosity *= gradSigmaEM;
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
        
        saturation = dielectricPermittiviyDeSaturation * gradEpsilonEM;
        if (model.getParameterisation() == 2) {
            // Based on Archie equation
            conductivityDeSaturation = this->getElectricConductivityDeSaturation(model); 
            conductivityDeSaturation *= gradSigmaEM;
            saturation += conductivityDeSaturation; 
        }           
    } else {
        this->initParameterisation(saturation, ctx, dist, 0.0);
    }  
    
    if (workflow.getInvertForReflectivity()) {
        reflectivity = -gradEpsilonEM;
    } else {
        this->initParameterisation(reflectivity, ctx, dist, 0.0);
    }  
}

/*! \brief calculate the stabilizing functional of each model parameter */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPrioriEM, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Workflow::Workflow<ValueType> const &workflow)
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
}

/*! \brief Apply a median filter to filter the extreme value of the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::EMEM<ValueType>::applyMedianFilter(KITGPI::Configuration::Configuration config, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> sigmaEM_temp;
    scai::lama::DenseVector<ValueType> epsilonEM_temp;
    scai::lama::DenseVector<ValueType> porositytemp;
    scai::lama::DenseVector<ValueType> saturationtemp;
    scai::lama::DenseVector<ValueType> reflectivitytemp;
    
    sigmaEM_temp = this->getElectricConductivity();
    epsilonEM_temp = this->getDielectricPermittivity();
    porositytemp = this->getPorosity();
    saturationtemp = this->getSaturation();
    reflectivitytemp = this->getReflectivity();
    
    scai::IndexType NZ = config.get<IndexType>("NZ");
    if (NZ == 1) {
        scai::IndexType NY = config.get<IndexType>("NY");
        scai::IndexType NX = porositytemp.size() / NY;
        scai::IndexType spatialLength = config.get<IndexType>("spatialFDorder");
        
        KITGPI::Common::applyMedianFilterTo2DVector(sigmaEM_temp, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(epsilonEM_temp, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(porositytemp, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(saturationtemp, NX, NY, spatialLength);
        KITGPI::Common::applyMedianFilterTo2DVector(reflectivitytemp, NX, NY, spatialLength);
        
        this->setElectricConductivity(sigmaEM_temp);    
        this->setDielectricPermittivity(epsilonEM_temp);
        this->setPorosity(porositytemp);    
        this->setSaturation(saturationtemp);
        this->setReflectivity(reflectivitytemp);
    }
}

template class KITGPI::Gradient::EMEM<float>;
template class KITGPI::Gradient::EMEM<double>;
