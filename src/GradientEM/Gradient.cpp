#include "Gradient.hpp"

using namespace scai;
using namespace KITGPI;

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Gradient::GradientEM<ValueType>::getRelaxationFrequency() const
{
    return (relaxationFrequency);
}

/*! \brief Set relaxation frequency
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setRelaxationFrequency(ValueType const setRelaxationFrequency)
{
    relaxationFrequency = setRelaxationFrequency;
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Gradient::GradientEM<ValueType>::getNumRelaxationMechanisms() const
{
    return (numRelaxationMechanisms);
}

/*! \brief Set number of Relaxation Mechanisms
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setNumRelaxationMechanisms(IndexType const setNumRelaxationMechanisms)
{
    numRelaxationMechanisms = setNumRelaxationMechanisms;
}

/*! \brief Init a single parameter by a constant value
 *
 \param vector Singel parameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param value Value which will be used to initialize the single parameter to a homogenoeus parameter vector
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value)
{
    allocateParameterisation(vector, ctx, dist);

    vector.setScalar(value);
}

/*! \brief Init a single parameter by reading a parameter vector from an external file
 *
 *  Reads a single parameter from an external mtx file.
 \param vector Singel parameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param filename Location of external file which will be read in
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{

    allocateParameterisation(vector, ctx, dist);

    readParameterisation(vector, filename, fileFormat);

    vector.redistribute(dist);
}

/*! \brief Write singe parameter to an external file
 *
 *  Write a single parameter to an external file block.
 \param vector Single parameter which will be written to filename
 \param filename Name of file in which parameter will be written
 \param fileFormat format for output file
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::writeParameterisation(scai::lama::Vector<ValueType> const &vector, std::string filename, IndexType fileFormat) const
{    
    IO::writeVector(vector, filename, fileFormat);
};

/*! \brief Read a parameter from file
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::readParameterisation(scai::lama::Vector<ValueType> &vector, std::string filename, IndexType fileFormat)
{    
    IO::readVector(vector, filename, fileFormat);
};

/*! \brief Allocate a single parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::allocateParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
}

/*! \brief If stream configuration is used, calculate a weighting vector to balance gradient per shot
 \param modelPerShot model per shot
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 \param uniqueShotInds unique shot indexes 
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::calcWeightingVector(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, std::vector<scai::IndexType> uniqueShotInds)
{
    auto dist = modelPerShot.getDielectricPermittivity().getDistributionPtr();
    auto distBig = dielectricPermittivity.getDistributionPtr();
    scai::hmemo::ContextPtr ctx = dielectricPermittivity.getContextPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix;
    scai::lama::CSRSparseMatrix<ValueType> recoverMatrix;
    shrinkMatrix.allocate(dist, distBig);
    recoverMatrix.allocate(distBig, dist);
    shrinkMatrix.setContextPtr(ctx);
    recoverMatrix.setContextPtr(ctx);
    scai::lama::DenseVector<ValueType> weightingVectorPerShot(dist, 1.0);
    scai::lama::DenseVector<ValueType> weightingVectorTemp(distBig, 0.0);
    IndexType numshots = uniqueShotInds.size();
    for (IndexType shotInd = 0; shotInd < numshots; shotInd++) {        
        shrinkMatrix = modelPerShot.getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(uniqueShotInds[shotInd]));
        recoverMatrix.assignTranspose(shrinkMatrix);
        weightingVectorTemp += recoverMatrix * weightingVectorPerShot;
    }
    weightingVector.allocate(distBig);
    weightingVector = 1.0 / weightingVectorTemp; // the weighting of overlapping area
    Common::replaceInvalid<ValueType>(weightingVector, 0.0);
    
    IO::writeVector(weightingVector, "gradients/weightingVectorEM", 1);
}

/*! \brief get weightingVector */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientEM<ValueType>::getWeightingVector()
{
    return weightingVector;
}
    
/*! \brief initialize an inner workflow */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setInvertForParameters(std::vector<bool> setInvertForParameters)
{
    workflowInner.setInvertForParameters(setInvertForParameters);
}
    
/*! \brief get an inner workflow */
template <typename ValueType>
std::vector<bool> KITGPI::Gradient::GradientEM<ValueType>::getInvertForParameters()
{
    return workflowInner.getInvertForParameters();
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
void KITGPI::Gradient::GradientEM<ValueType>::gradientParameterisation(scai::lama::DenseVector<ValueType> &vecGradientParameter, scai::lama::DenseVector<ValueType> vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation)
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
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getElectricConductivityDePorosity(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model)
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
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getDielectricPermittiviyDePorosity(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model)
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
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getElectricConductivityDeSaturation(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model)
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
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getDielectricPermittiviyDeSaturation(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model)
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

template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration config, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> electricConductivitytemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivitytemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempX;
    scai::lama::DenseVector<ValueType> tempY;
    ValueType modelMax;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(config, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(config, "NY");
    scai::IndexType spatialFDorder = config.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = config.get<IndexType>("exchangeStrategy");
    ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum();
    ValueType const ElectricConductivityReference = model.getElectricConductivityReference();
        
    /* Get references to required derivatives matrixes */
    auto const &DxfEM = derivativesEM.getDxf();
    auto const &DyfEM = derivativesEM.getDyf();  
          
    electricConductivitytemp = model.getElectricConductivity();  
    dielectricPermittivitytemp = model.getDielectricPermittivity();             
                      
    scai::hmemo::ContextPtr ctx = dielectricPermittivitytemp.getContextPtr();
    scai::dmemo::DistributionPtr dist = dielectricPermittivitytemp.getDistributionPtr();   
      
    this->applyParameterisation(electricConductivitytemp, ElectricConductivityReference, model.getParameterisation());       
    this->applyParameterisation(dielectricPermittivitytemp, DielectricPermittivityVacuum, model.getParameterisation());       
               
    if (workflow.getInvertForEpsilonEM()) {
        tempX = DxfEM * dielectricPermittivitytemp;    
        tempY = DyfEM * dielectricPermittivitytemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder);  
        
        modelMax = tempX.maxNorm();
        if (modelMax != 0)
            tempX *= 1 / modelMax;  
        modelMax = tempY.maxNorm();
        if (modelMax != 0)
            tempY *= 1 / modelMax;                       
    } else {
        this->initParameterisation(tempX, ctx, dist, 0.0);
        this->initParameterisation(tempY, ctx, dist, 0.0);
    }    
    
    modelDerivativeXtemp = tempX;
    modelDerivativeYtemp = tempY;   
    
    if (exchangeStrategy == 2) {
        if (workflow.getInvertForSigmaEM()) {        
            tempX = DxfEM * electricConductivitytemp;    
            tempY = DyfEM * electricConductivitytemp;
            
            KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder);
            KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder);   
            
            modelMax = tempX.maxNorm();
            if (modelMax != 0)
                tempX *= 1 / modelMax;  
            modelMax = tempY.maxNorm();
            if (modelMax != 0)
                tempY *= 1 / modelMax;                                       
        } else {
            this->initParameterisation(tempX, ctx, dist, 0.0);
            this->initParameterisation(tempY, ctx, dist, 0.0);
        }   
        
        modelDerivativeXtemp += tempX;
        modelDerivativeYtemp += tempY; 
    }
    
    dataMisfitEM.setModelDerivativeX(modelDerivativeXtemp);
    dataMisfitEM.setModelDerivativeY(modelDerivativeYtemp);
}

template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration config, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> electricConductivitytemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivitytemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> temp;
    scai::lama::DenseVector<ValueType> tempX;
    scai::lama::DenseVector<ValueType> tempY;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(config, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(config, "NY");
    scai::IndexType spatialFDorder = config.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = config.get<IndexType>("exchangeStrategy");
    ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum();
    ValueType const ElectricConductivityReference = model.getElectricConductivityReference();
    
    /* Get references to required derivatives matrixes */
    auto const &DxfEM = derivativesEM.getDxf();
    auto const &DyfEM = derivativesEM.getDyf();  
    
    electricConductivitytemp = model.getElectricConductivity();
    dielectricPermittivitytemp = model.getDielectricPermittivity();  
        
    this->applyParameterisation(electricConductivitytemp, ElectricConductivityReference, model.getParameterisation());       
    this->applyParameterisation(dielectricPermittivitytemp, DielectricPermittivityVacuum, model.getParameterisation());       
    
    scai::hmemo::ContextPtr ctx = dielectricPermittivitytemp.getContextPtr();
    scai::dmemo::DistributionPtr dist = dielectricPermittivitytemp.getDistributionPtr();        
                              
    modelDerivativeXtemp = dataMisfit.getModelDerivativeX();  
    modelDerivativeYtemp = dataMisfit.getModelDerivativeY();   
    modelDerivativeXtemp = modelTaper2DJoint.applyGradientTransformToEM(modelDerivativeXtemp); 
    modelDerivativeYtemp = modelTaper2DJoint.applyGradientTransformToEM(modelDerivativeYtemp); 
                 
    if (workflow.getInvertForEpsilonEM()) {    
        tempX = DxfEM * dielectricPermittivitytemp;    
        tempY = DyfEM * dielectricPermittivitytemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder); 
        
        if (exchangeStrategy != 0) {
            // cross-gradient of vs and modelDerivative   
            dielectricPermittivitytemp = modelDerivativeYtemp * tempX;   
            temp = modelDerivativeXtemp * tempY;   
            dielectricPermittivitytemp -= temp;  
            
            dielectricPermittivity = dielectricPermittivitytemp;   
        }
    } else {
        this->initParameterisation(dielectricPermittivity, ctx, dist, 0.0);
        this->initParameterisation(tempX, ctx, dist, 0.0);
        this->initParameterisation(tempY, ctx, dist, 0.0);
    }   
    
    modelDerivativeXtemp = tempX;
    modelDerivativeYtemp = tempY; 
            
    if (workflow.getInvertForSigmaEM()) {
        tempX = DxfEM * electricConductivitytemp;    
        tempY = DyfEM * electricConductivitytemp;          
        
        // cross-gradient of electricConductivity and modelDerivative   
        electricConductivitytemp = modelDerivativeYtemp * tempX;   
        temp = modelDerivativeXtemp * tempY;   
        electricConductivitytemp -= temp;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(electricConductivitytemp, NX, NY, spatialFDorder);        
        electricConductivity = electricConductivitytemp;                
    } else {
        this->initParameterisation(electricConductivity, ctx, dist, 0.0);
    }       
}

template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration config, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow)
{
    scai::lama::DenseVector<ValueType> electricConductivitytemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivitytemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempX;
    scai::lama::DenseVector<ValueType> tempY;
    scai::lama::DenseVector<ValueType> temp;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(config, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(config, "NY");
    scai::IndexType spatialFDorder = config.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = config.get<IndexType>("exchangeStrategy");
    ValueType const DielectricPermittivityVacuum = model.getDielectricPermittivityVacuum();
        
    dielectricPermittivitytemp = model.getDielectricPermittivity();
    this->applyParameterisation(dielectricPermittivitytemp, DielectricPermittivityVacuum, model.getParameterisation());      
               
    modelDerivativeXtemp = dataMisfit.getModelDerivativeX();  
    modelDerivativeYtemp = dataMisfit.getModelDerivativeY();
                      
    scai::hmemo::ContextPtr ctx = dielectricPermittivitytemp.getContextPtr();
    scai::dmemo::DistributionPtr dist = dielectricPermittivitytemp.getDistributionPtr();  
        
    /* Get references to required derivatives matrixes */
    auto const &DxfEM = derivativesEM.getDxf();
    auto const &DyfEM = derivativesEM.getDyf();  
    
    modelDerivativeXtemp = modelTaper2DJoint.applyGradientTransformToEM(modelDerivativeXtemp); 
    modelDerivativeYtemp = modelTaper2DJoint.applyGradientTransformToEM(modelDerivativeYtemp);   
            
    temp = modelDerivativeXtemp - modelDerivativeYtemp; 
                 
    if (workflow.getInvertForEpsilonEM()) {
        tempX = DxfEM * dielectricPermittivitytemp;    
        tempY = DyfEM * dielectricPermittivitytemp;   
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder);
        
        if (exchangeStrategy != 0) {
            // derivative of cross-gradient with respect to dielectricPermittivity 
            dielectricPermittivitytemp = this->getDielectricPermittivity();
            dielectricPermittivitytemp *= temp;  
            
            dielectricPermittivity = dielectricPermittivitytemp;          
        } else {
            this->initParameterisation(dielectricPermittivity, ctx, dist, 0.0);
        }              
    } else {
        this->initParameterisation(dielectricPermittivity, ctx, dist, 0.0);
        this->initParameterisation(tempX, ctx, dist, 0.0);
        this->initParameterisation(tempY, ctx, dist, 0.0);
    }        
    
    temp = tempX - tempY; 
        
    if (workflow.getInvertForSigmaEM()) {      
        // derivative of cross-gradient with respect to electricConductivity   
        electricConductivitytemp = this->getElectricConductivity();
        electricConductivitytemp *= temp;  
                      
        electricConductivity = electricConductivitytemp;
    } else {
        this->initParameterisation(electricConductivity, ctx, dist, 0.0);
    }     
}

/*! \brief calculate the gradient of stabilizing functional of each model parameter */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientEM<ValueType>::calcStabilizingFunctionalGradientPerModel(scai::lama::DenseVector<ValueType> modelResidualVec, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM)
{ 
    scai::lama::DenseVector<ValueType> regularizedGradient; 
    scai::lama::DenseVector<ValueType> We; 
    ValueType tempValue;
    tempValue = config.get<ValueType>("focusingParameter");
    tempValue *= modelResidualVec.maxNorm();
    tempValue = pow(tempValue, 2);  
    We = scai::lama::pow(modelResidualVec, 2);
    We += tempValue;
    tempValue = dataMisfitEM.calcStablizingFunctionalPerModel(modelResidualVec, config.get<ValueType>("focusingParameter"), config.get<ValueType>("stablizingFunctionalType"));
    We = tempValue / We;
    We = scai::lama::sqrt(We);   
    
    regularizedGradient = We * We; 
    regularizedGradient *= modelResidualVec;
    
    return regularizedGradient;
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

/*! \brief Get const reference to porosity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getPorosity()
{
    return (porosity);
}

/*! \brief Get const reference to porosity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getPorosity() const
{
    return (porosity);
}

/*! \brief Set porosity
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setPorosity(scai::lama::Vector<ValueType> const &setPorosity)
{
    porosity = setPorosity;
}

/*! \brief Get const reference to saturation
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getSaturation()
{
    return (saturation);
}

/*! \brief Get const reference to saturation
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getSaturation() const
{
    return (saturation);
}

/*! \brief Set saturation
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setSaturation(scai::lama::Vector<ValueType> const &setSaturation)
{
    saturation = setSaturation;
}

/*! \brief Get const reference to reflectivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getReflectivity()
{
    return (reflectivity);
}

/*! \brief Get const reference to reflectivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getReflectivity() const
{
    return (reflectivity);
}

/*! \brief Set reflectivity
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setReflectivity(scai::lama::Vector<ValueType> const &setReflectivity)
{
    reflectivity = setReflectivity;
}

/*! \brief Get parameter normalizeGradient
 */
template <typename ValueType>
bool KITGPI::Gradient::GradientEM<ValueType>::getNormalizeGradient() const
{
    return (normalizeGradient);
}

/*! \brief Set normalizeGradient parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setNormalizeGradient(bool const &gradNorm)
{
    normalizeGradient = gradNorm;
}

/*! \brief Overloading = Operation
 *
 \param rhs Gradient which is copied.
 */
template <typename ValueType>
KITGPI::Gradient::GradientEM<ValueType> &KITGPI::Gradient::GradientEM<ValueType>::operator=(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    assign(rhs);
    return *this;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::GradientEM<ValueType> &KITGPI::Gradient::GradientEM<ValueType>::operator-=(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    minusAssign(rhs);
    return *this;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::GradientEM<ValueType> &KITGPI::Gradient::GradientEM<ValueType>::operator+=(KITGPI::Gradient::GradientEM<ValueType> const &rhs)
{
    plusAssign(rhs);
    return *this;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::GradientEM<ValueType> &KITGPI::Gradient::GradientEM<ValueType>::operator*=(ValueType const &rhs)
{
    timesAssign(rhs);
    return *this;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::GradientEM<ValueType> &KITGPI::Gradient::GradientEM<ValueType>::operator*=(scai::lama::Vector<ValueType> const &rhs)
{
    timesAssign(rhs);
    return *this;
}

template class KITGPI::Gradient::GradientEM<float>;
template class KITGPI::Gradient::GradientEM<double>;
