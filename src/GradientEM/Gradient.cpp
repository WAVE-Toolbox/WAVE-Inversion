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
 \param distEM Distribution
 \param value Value which will be used to initialize the single parameter to a homogenoeus parameter vector
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, ValueType value)
{
    allocateParameterisation(vector, ctx, distEM);

    vector.setScalar(value);
}

/*! \brief Init a single parameter by reading a parameter vector from an external file
 *
 *  Reads a single parameter from an external mtx file.
 \param vector Singel parameter which will be initialized
 \param ctx Context
 \param distEM Distribution
 \param filename Location of external file which will be read in
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, std::string filename, IndexType fileFormat)
{

    allocateParameterisation(vector, ctx, distEM);

    readParameterisation(vector, filename, fileFormat);

    vector.redistribute(distEM);
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
void KITGPI::Gradient::GradientEM<ValueType>::allocateParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM)
{
    vector.setContextPtr(ctx);
    vector.allocate(distEM);
};

/*! \brief Get matrix that multiplies with model matrices to get a pershot
 \param dist Distribution of the pershot
 \param distBig Distribution of the big model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate coordinate where to cut the pershot
 */
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> KITGPI::Gradient::GradientEM<ValueType>::getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
{
    SparseFormat shrinkMatrix; //!< Shrink Multiplication matrix
    shrinkMatrix.allocate(dist, distBig);
      
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);
    lama::MatrixAssembly<ValueType> assembly;
    IndexType indexBig = 0;
 
    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) { // loop over all indices
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex); // get submodel coordinate from singleIndex
        coordinate.x += cutCoordinate.x; // offset depends on shot number
        coordinate.y += cutCoordinate.y;
        coordinate.z += cutCoordinate.z;
        indexBig = modelCoordinatesBig.coordinate2index(coordinate);
        assembly.push(ownedIndex, indexBig, 1.0);
    }
    shrinkMatrix.fillFromAssembly(assembly);
    return shrinkMatrix;
}

/*! \brief Get erase-matrix that erases the old values in the big model
 \param dist Distribution of the pershot
 \param distBig Distribution of the big model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate coordinate where to cut the pershot
 */
template <typename ValueType>
scai::lama::SparseVector<ValueType> KITGPI::Gradient::GradientEM<ValueType>::getEraseVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth)
{
    scai::lama::SparseVector<ValueType> eraseVector(distBig, 1.0); //!< Shrink Multiplication matrix
      
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);
    lama::VectorAssembly<ValueType> assembly;
    IndexType indexBig = 0;
 
    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) { // loop over all indices
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex); // get submodel coordinate from singleIndex
        coordinate.x += cutCoordinate.x; // offset depends on shot number
        coordinate.y += cutCoordinate.y;
        coordinate.z += cutCoordinate.z;
        indexBig = modelCoordinatesBig.coordinate2index(coordinate);
        assembly.push(indexBig, 0.0);
    }
    eraseVector.fillFromAssembly(assembly);
    
    // damp the boundary boarders
    for (IndexType y = 0; y < modelCoordinatesBig.getNY(); y++) {
        for (IndexType i = 0; i < boundaryWidth; i++) {
            ValueType tmp = (ValueType)1.0 - (i + 1) / (ValueType)boundaryWidth;
            eraseVector[modelCoordinatesBig.coordinate2index(cutCoordinate.x+i, y, 0)] = tmp;
            eraseVector[modelCoordinatesBig.coordinate2index(cutCoordinate.x+modelCoordinates.getNX()-1-i, y, 0)] = tmp;
        }
    }
    return eraseVector;
}

/*! \brief initialize an inner workflow */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setInvertParameterSingle(std::vector<bool> setInvertParameterSingle)
{
    workflowInner.setInvertParameters(setInvertParameterSingle);
}
    
/*! \brief get an inner workflow */
template <typename ValueType>
std::vector<bool> KITGPI::Gradient::GradientEM<ValueType>::getInvertParameterSingle()
{
    return workflowInner.getInvertParameters();
}

/*! \brief apply parameterisation of modelEM parameters */
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

/*! \brief execute parameterisation of modelEM parameters */
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

/*! \brief delete parameterisation of modelEM parameters */
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
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getConductivityDePorosity(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM)
{
    scai::lama::DenseVector<ValueType> conductivityDePorosity;
    scai::lama::DenseVector<ValueType> temp;  

    // \pdv{\sigma_e}{\phi} = \frac{m}{a} \sigma_{ew} \phi^{m-1} S_w^{n}
    ValueType aArchietemp;
    ValueType mArchietemp;
    ValueType nArchietemp;
    aArchietemp = modelEM.getArchie_a();
    mArchietemp = modelEM.getArchie_m();
    nArchietemp = modelEM.getArchie_n();
    // Based on Archie equation
    conductivityDePorosity = modelEM.getConductivityEMWater();   
    conductivityDePorosity *= mArchietemp / aArchietemp;
    temp = scai::lama::pow(modelEM.getPorosity(), mArchietemp - 1); 
    Common::replaceInvalid<ValueType>(temp, 0.0);
    conductivityDePorosity *= temp;
    temp = scai::lama::pow(modelEM.getSaturation(), nArchietemp); 
    conductivityDePorosity *= temp; 
    
    return conductivityDePorosity;
}

/*! \brief calculate the derivative of dielectricPermittiviy with respect to porosity */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getDielectricPermittiviyDePorosity(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM)
{
    scai::lama::DenseVector<ValueType> dielectricPermittiviyDePorosity;
    scai::lama::DenseVector<ValueType> temp; 
    ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum(); 
    ValueType const RelativeDielectricPermittivityWater = modelEM.getRelativeDielectricPermittivityWater(); 
    ValueType const RelativeDielectricPermittivityVacuum = modelEM.getRelativeDielectricPermittivityVacuum(); 
    
    // \pdv{\varepsilon_{e}}{\phi} = 2 \left[ - \sqrt{ \varepsilon_{emar} } + S_w \sqrt{ \varepsilon_{ewr} } + \left(1-S_w\right) \sqrt{ \varepsilon_{ear}} \right] \sqrt{ \varepsilon_{e} \varepsilon_{e0}}
    // Based on complex refractive index model (CRIM)    
    dielectricPermittiviyDePorosity = scai::lama::sqrt(modelEM.getRelativeDieletricPeimittivityRockMatrix()); 
    temp = modelEM.getSaturation();  
    temp *= sqrt(RelativeDielectricPermittivityWater);
    dielectricPermittiviyDePorosity = temp - dielectricPermittiviyDePorosity;
    temp = 1 - modelEM.getSaturation();  
    temp *= sqrt(RelativeDielectricPermittivityVacuum);
    dielectricPermittiviyDePorosity += temp;    
    temp = modelEM.getDielectricPermittivityEM();   
    temp *= DielectricPermittivityVacuum;
    temp = scai::lama::sqrt(temp); 
    dielectricPermittiviyDePorosity *= temp;   
    dielectricPermittiviyDePorosity *= 2; 
    
    return dielectricPermittiviyDePorosity;
}

/*! \brief calculate the derivative of conductivity with respect to saturation */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getConductivityDeSaturation(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM)
{
    scai::lama::DenseVector<ValueType> conductivityDeSaturation;
    scai::lama::DenseVector<ValueType> temp;  
    
    // \pdv{\sigma_e}{S_w} = \frac{n}{a} \sigma_{ew} \phi^m S_w^{n-1}
    ValueType aArchietemp;
    ValueType mArchietemp;
    ValueType nArchietemp;
    aArchietemp = modelEM.getArchie_a();
    mArchietemp = modelEM.getArchie_m();
    nArchietemp = modelEM.getArchie_n();
    // Based on Archie equation
    conductivityDeSaturation = modelEM.getConductivityEMWater();   
    conductivityDeSaturation *= nArchietemp / aArchietemp;
    temp = scai::lama::pow(modelEM.getPorosity(), mArchietemp); 
    conductivityDeSaturation *= temp;
    temp = scai::lama::pow(modelEM.getSaturation(), nArchietemp - 1); 
    Common::replaceInvalid<ValueType>(temp, 0.0);
    conductivityDeSaturation *= temp; 
    
    return conductivityDeSaturation;
}

/*! \brief calculate the derivative of conductivity with respect to saturation */
template <typename ValueType> scai::lama::DenseVector<ValueType>  KITGPI::Gradient::GradientEM<ValueType>::getDielectricPermittiviyDeSaturation(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM)
{
    scai::lama::DenseVector<ValueType> dielectricPermittiviyDeSaturation;
    scai::lama::DenseVector<ValueType> temp;  
    ValueType tempValue;
    ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
    ValueType const RelativeDielectricPermittivityWater = modelEM.getRelativeDielectricPermittivityWater(); 
    ValueType const RelativeDielectricPermittivityVacuum = modelEM.getRelativeDielectricPermittivityVacuum(); 
        
    // \pdv{\varepsilon_{e}}{S_w} = 2  \phi \left( \sqrt{ \varepsilon_{ewr} } - \sqrt{ \varepsilon_{ear}} \right) \sqrt{ \varepsilon_{e} \varepsilon_{e0}}
    // Based on complex refractive index model (CRIM)    
    temp = modelEM.getDielectricPermittivityEM();   
    temp *= DielectricPermittivityVacuum;
    dielectricPermittiviyDeSaturation = scai::lama::sqrt(temp);         
    temp = modelEM.getPorosity();  
    tempValue = sqrt(RelativeDielectricPermittivityWater) - sqrt(RelativeDielectricPermittivityVacuum);
    temp *= tempValue;
    dielectricPermittiviyDeSaturation *= temp;   
    dielectricPermittiviyDeSaturation *= 2;  
    
    return dielectricPermittiviyDeSaturation;
}

template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    scai::lama::DenseVector<ValueType> conductivityEMtemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEMtemp;
    scai::lama::DenseVector<ValueType> modelEMDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelEMDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempX;
    scai::lama::DenseVector<ValueType> tempY;
    ValueType modelMax;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialFDorder = configEM.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = configEM.get<IndexType>("exchangeStrategy");
    ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
    ValueType const ConductivityReference = modelEM.getConductivityReference();
        
    /* Get references to required derivatives matrixes */
    auto const &DxfEM = derivativesEM.getDxf();
    auto const &DyfEM = derivativesEM.getDyf();  
          
    conductivityEMtemp = modelEM.getConductivityEM();  
    dielectricPermittivityEMtemp = modelEM.getDielectricPermittivityEM();             
                      
    scai::hmemo::ContextPtr ctx = dielectricPermittivityEMtemp.getContextPtr();
    scai::dmemo::DistributionPtr distEM = dielectricPermittivityEMtemp.getDistributionPtr();   
      
    this->applyParameterisation(conductivityEMtemp, ConductivityReference, modelEM.getParameterisation());       
    this->applyParameterisation(dielectricPermittivityEMtemp, DielectricPermittivityVacuum, modelEM.getParameterisation());       
               
    if (workflowEM.getInvertForEpsilonEM()) {
        tempX = DxfEM * dielectricPermittivityEMtemp;    
        tempY = DyfEM * dielectricPermittivityEMtemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder);
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder);  
        
        modelMax = tempX.maxNorm();
        if (modelMax != 0)
            tempX *= 1 / modelMax;  
        modelMax = tempY.maxNorm();
        if (modelMax != 0)
            tempY *= 1 / modelMax;                       
    } else {
        this->initParameterisation(tempX, ctx, distEM, 0.0);
        this->initParameterisation(tempY, ctx, distEM, 0.0);
    }    
    
    modelEMDerivativeXtemp = tempX;
    modelEMDerivativeYtemp = tempY;   
    
    if (exchangeStrategy == 2) {
        if (workflowEM.getInvertForSigmaEM()) {        
            tempX = DxfEM * conductivityEMtemp;    
            tempY = DyfEM * conductivityEMtemp;
            
            KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder);
            KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder);   
            
            modelMax = tempX.maxNorm();
            if (modelMax != 0)
                tempX *= 1 / modelMax;  
            modelMax = tempY.maxNorm();
            if (modelMax != 0)
                tempY *= 1 / modelMax;                                       
        } else {
            this->initParameterisation(tempX, ctx, distEM, 0.0);
            this->initParameterisation(tempY, ctx, distEM, 0.0);
        }   
        
        modelEMDerivativeXtemp += tempX;
        modelEMDerivativeYtemp += tempY; 
    }
    
    dataMisfitEM.setModelDerivativeX(modelEMDerivativeXtemp);
    dataMisfitEM.setModelDerivativeY(modelEMDerivativeYtemp);
}

template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    scai::lama::DenseVector<ValueType> conductivityEMtemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEMtemp;
    scai::lama::DenseVector<ValueType> modelEMDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelEMDerivativeYtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> temp;
    scai::lama::DenseVector<ValueType> tempX;
    scai::lama::DenseVector<ValueType> tempY;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialFDorder = configEM.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = configEM.get<IndexType>("exchangeStrategy");
    ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
    ValueType const ConductivityReference = modelEM.getConductivityReference();
    
    /* Get references to required derivatives matrixes */
    auto const &DxfEM = derivativesEM.getDxf();
    auto const &DyfEM = derivativesEM.getDyf();  
    
    conductivityEMtemp = modelEM.getConductivityEM();
    dielectricPermittivityEMtemp = modelEM.getDielectricPermittivityEM();  
        
    this->applyParameterisation(conductivityEMtemp, ConductivityReference, modelEM.getParameterisation());       
    this->applyParameterisation(dielectricPermittivityEMtemp, DielectricPermittivityVacuum, modelEM.getParameterisation());       
    
    scai::hmemo::ContextPtr ctx = dielectricPermittivityEMtemp.getContextPtr();
    scai::dmemo::DistributionPtr distEM = dielectricPermittivityEMtemp.getDistributionPtr();        
                              
    modelDerivativeXtemp = dataMisfit.getModelDerivativeX();  
    modelDerivativeYtemp = dataMisfit.getModelDerivativeY();   
    modelDerivativeXtemp = modelTaper2DJoint.applyGradientTransformToEM(modelDerivativeXtemp); 
    modelDerivativeYtemp = modelTaper2DJoint.applyGradientTransformToEM(modelDerivativeYtemp); 
                 
    if (workflowEM.getInvertForEpsilonEM()) {    
        tempX = DxfEM * dielectricPermittivityEMtemp;    
        tempY = DyfEM * dielectricPermittivityEMtemp;
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder); 
        
        if (exchangeStrategy != 0) {
            // cross-gradient of vs and modelDerivative   
            dielectricPermittivityEMtemp = modelDerivativeYtemp * tempX;   
            temp = modelDerivativeXtemp * tempY;   
            dielectricPermittivityEMtemp -= temp;  
            
            dielectricPermittivityEM = dielectricPermittivityEMtemp;   
        }
    } else {
        this->initParameterisation(dielectricPermittivityEM, ctx, distEM, 0.0);
        this->initParameterisation(tempX, ctx, distEM, 0.0);
        this->initParameterisation(tempY, ctx, distEM, 0.0);
    }   
    
    modelEMDerivativeXtemp = tempX;
    modelEMDerivativeYtemp = tempY; 
            
    if (workflowEM.getInvertForSigmaEM()) {
        tempX = DxfEM * conductivityEMtemp;    
        tempY = DyfEM * conductivityEMtemp;          
        
        // cross-gradient of conductivityEM and modelDerivative   
        conductivityEMtemp = modelEMDerivativeYtemp * tempX;   
        temp = modelEMDerivativeXtemp * tempY;   
        conductivityEMtemp -= temp;  
        
        KITGPI::Common::applyMedianFilterTo2DVector(conductivityEMtemp, NX, NY, spatialFDorder);        
        conductivityEM = conductivityEMtemp;                
    } else {
        this->initParameterisation(conductivityEM, ctx, distEM, 0.0);
    }       
}

template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    scai::lama::DenseVector<ValueType> conductivityEMtemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEMtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeXtemp;
    scai::lama::DenseVector<ValueType> modelDerivativeYtemp;
    scai::lama::DenseVector<ValueType> tempX;
    scai::lama::DenseVector<ValueType> tempY;
    scai::lama::DenseVector<ValueType> temp;
    scai::IndexType NX = Common::getFromStreamFile<IndexType>(configEM, "NX");
    scai::IndexType NY = Common::getFromStreamFile<IndexType>(configEM, "NY");
    scai::IndexType spatialFDorder = configEM.get<IndexType>("spatialFDorder");
    scai::IndexType exchangeStrategy = configEM.get<IndexType>("exchangeStrategy");
    ValueType const DielectricPermittivityVacuum = modelEM.getDielectricPermittivityVacuum();
        
    dielectricPermittivityEMtemp = modelEM.getDielectricPermittivityEM();
    this->applyParameterisation(dielectricPermittivityEMtemp, DielectricPermittivityVacuum, modelEM.getParameterisation());      
               
    modelDerivativeXtemp = dataMisfit.getModelDerivativeX();  
    modelDerivativeYtemp = dataMisfit.getModelDerivativeY();
                      
    scai::hmemo::ContextPtr ctx = dielectricPermittivityEMtemp.getContextPtr();
    scai::dmemo::DistributionPtr distEM = dielectricPermittivityEMtemp.getDistributionPtr();  
        
    /* Get references to required derivatives matrixes */
    auto const &DxfEM = derivativesEM.getDxf();
    auto const &DyfEM = derivativesEM.getDyf();  
    
    modelDerivativeXtemp = modelTaper2DJoint.applyGradientTransformToEM(modelDerivativeXtemp); 
    modelDerivativeYtemp = modelTaper2DJoint.applyGradientTransformToEM(modelDerivativeYtemp);   
            
    temp = modelDerivativeXtemp - modelDerivativeYtemp; 
                 
    if (workflowEM.getInvertForEpsilonEM()) {
        tempX = DxfEM * dielectricPermittivityEMtemp;    
        tempY = DyfEM * dielectricPermittivityEMtemp;   
        
        KITGPI::Common::applyMedianFilterTo2DVector(tempX, NX, NY, spatialFDorder); 
        KITGPI::Common::applyMedianFilterTo2DVector(tempY, NX, NY, spatialFDorder);
        
        if (exchangeStrategy != 0) {
            // derivative of cross-gradient with respect to dielectricPermittivityEM 
            dielectricPermittivityEMtemp = this->getDielectricPermittivityEM();
            dielectricPermittivityEMtemp *= temp;  
            
            dielectricPermittivityEM = dielectricPermittivityEMtemp;          
        } else {
            this->initParameterisation(dielectricPermittivityEM, ctx, distEM, 0.0);
        }              
    } else {
        this->initParameterisation(dielectricPermittivityEM, ctx, distEM, 0.0);
        this->initParameterisation(tempX, ctx, distEM, 0.0);
        this->initParameterisation(tempY, ctx, distEM, 0.0);
    }        
    
    temp = tempX - tempY; 
        
    if (workflowEM.getInvertForSigmaEM()) {      
        // derivative of cross-gradient with respect to conductivityEM   
        conductivityEMtemp = this->getConductivityEM();
        conductivityEMtemp *= temp;  
                      
        conductivityEM = conductivityEMtemp;
    } else {
        this->initParameterisation(conductivityEM, ctx, distEM, 0.0);
    }     
}

/*! \brief calculate the gradient of stabilizing functional of each model parameter */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::GradientEM<ValueType>::calcStabilizingFunctionalGradientPerModel(scai::lama::DenseVector<ValueType> modelResidualVec, KITGPI::Configuration::Configuration configEM, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM)
{ 
    scai::lama::DenseVector<ValueType> regularizedGradient; 
    scai::lama::DenseVector<ValueType> We; 
    ValueType tempValue;
    tempValue = configEM.get<ValueType>("focusingParameter");
    tempValue *= modelResidualVec.maxNorm();
    tempValue = pow(tempValue, 2);  
    We = scai::lama::pow(modelResidualVec, 2);
    We += tempValue;
    tempValue = dataMisfitEM.calcStablizingFunctionalPerModel(modelResidualVec, configEM.get<ValueType>("focusingParameter"), configEM.get<ValueType>("stablizingFunctionalType"));
    We = tempValue / We;
    We = scai::lama::sqrt(We);   
    
    regularizedGradient = We * We; 
    regularizedGradient *= modelResidualVec;
    
    return regularizedGradient;
}

/*! \brief Get const reference to conductivityEM
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getConductivityEM()
{
    return (conductivityEM);
}

/*! \brief Get const reference to conductivityEM
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getConductivityEM() const
{
    return (conductivityEM);
}

/*! \brief Set conductivityEM
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setConductivityEM(scai::lama::Vector<ValueType> const &setConductivityEM)
{
    conductivityEM = setConductivityEM;
}


/*! \brief Get const reference to dielectricPermittivityEM
 * 
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getDielectricPermittivityEM()
{
    return (dielectricPermittivityEM);
}

/*! \brief Get const reference to dielectricPermittivityEM
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getDielectricPermittivityEM() const
{
    return (dielectricPermittivityEM);
}

/*! \brief Set dielectricPermittivityEM
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setDielectricPermittivityEM(scai::lama::Vector<ValueType> const &setDielectricPermittivityEM)
{
    dielectricPermittivityEM = setDielectricPermittivityEM;
}

/*! \brief Get const reference to tauConductivityEM
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauConductivityEM()
{
    return (tauConductivityEM);
}

/*! \brief Get const reference to tauConductivityEM
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauConductivityEM() const
{
    return (tauConductivityEM);
}

/*! \brief Set tauConductivityEM parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setTauConductivityEM(scai::lama::Vector<ValueType> const &setTauConductivityEM)
{
    tauConductivityEM = setTauConductivityEM;
}

/*! \brief Get const reference to tauDielectricPermittivityEM
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauDielectricPermittivityEM()
{
    return (tauDielectricPermittivityEM);
}

/*! \brief Get const reference to tauDielectricPermittivityEM
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::GradientEM<ValueType>::getTauDielectricPermittivityEM() const
{
    return (tauDielectricPermittivityEM);
}

/*! \brief Set tauDielectricPermittivityEM parameter
 */
template <typename ValueType>
void KITGPI::Gradient::GradientEM<ValueType>::setTauDielectricPermittivityEM(scai::lama::Vector<ValueType> const &setTauDielectricPermittivityEM)
{
    tauDielectricPermittivityEM = setTauDielectricPermittivityEM;
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
