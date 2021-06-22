#include "Gradient.hpp"

using namespace scai;
using namespace KITGPI;

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Gradient::Gradient<ValueType>::getRelaxationFrequency() const
{
    return (relaxationFrequency);
}

/*! \brief Set relaxation frequency
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setRelaxationFrequency(ValueType const setRelaxationFrequency)
{
    relaxationFrequency = setRelaxationFrequency;
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Gradient::Gradient<ValueType>::getNumRelaxationMechanisms() const
{
    return (numRelaxationMechanisms);
}

/*! \brief Set number of Relaxation Mechanisms
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setNumRelaxationMechanisms(IndexType const setNumRelaxationMechanisms)
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
void KITGPI::Gradient::Gradient<ValueType>::initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value)
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
void KITGPI::Gradient::Gradient<ValueType>::initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
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
void KITGPI::Gradient::Gradient<ValueType>::writeParameterisation(scai::lama::Vector<ValueType> const &vector, std::string filename, IndexType fileFormat) const
{    
    IO::writeVector(vector, filename, fileFormat);
};

/*! \brief Read a parameter from file
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::readParameterisation(scai::lama::Vector<ValueType> &vector, std::string filename, IndexType fileFormat)
{    
    IO::readVector(vector, filename, fileFormat);
};

/*! \brief Allocate a single parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::allocateParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

/*! \brief Get matrix that multiplies with model matrices to get a pershot
 \param dist Distribution of the pershot
 \param distBig Distribution of the big model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate coordinate where to cut the pershot
 */
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> KITGPI::Gradient::Gradient<ValueType>::getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
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
scai::lama::SparseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::getEraseVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth)
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
void KITGPI::Gradient::Gradient<ValueType>::setInvertParameters(std::vector<bool> setInvertParameters)
{
    workflowInner.setInvertParameters(setInvertParameters);
}
    
/*! \brief get an inner workflow */
template <typename ValueType>
std::vector<bool> KITGPI::Gradient::Gradient<ValueType>::getInvertParameters()
{
    return workflowInner.getInvertParameters();
}

/*! \brief calculate the derivative of density with respect to porosity
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::getDensityDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
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
scai::lama::DenseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::getMu_satDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
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
scai::lama::DenseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::getK_satDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
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
scai::lama::DenseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::getBiotCoefficientDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
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
scai::lama::DenseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::getDensityDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
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
scai::lama::DenseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::getK_satDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
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

/*! \brief calculate the gradient of stabilizing functional of each model parameter */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::calcStabilizingFunctionalGradientPerModel(scai::lama::DenseVector<ValueType> modelResidualVec, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit)
{ 
    scai::lama::DenseVector<ValueType> regularizedGradient; 
    scai::lama::DenseVector<ValueType> We; 
    ValueType tempValue;
    tempValue = config.get<ValueType>("focusingParameter");
    tempValue *= modelResidualVec.maxNorm();
    tempValue = pow(tempValue, 2);  
    We = scai::lama::pow(modelResidualVec, 2);
    We += tempValue;
    tempValue = dataMisfit.calcStablizingFunctionalPerModel(modelResidualVec, config.get<ValueType>("focusingParameter"), config.get<ValueType>("stablizingFunctionalType"));
    We = tempValue / We;
    We = scai::lama::sqrt(We);   
    
    regularizedGradient = We * We; 
    regularizedGradient *= modelResidualVec;
    
    return regularizedGradient;
}

/*! \brief Get const reference to density parameter
 * 
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getDensity()
{

    return (density);
}

/*! \brief Get const reference to density parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getDensity() const
{
    return (density);
}

/*! \brief Set density  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setDensity(scai::lama::Vector<ValueType> const &setDensity)
{
    density = setDensity;
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getVelocityP()
{
    return (velocityP);
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getVelocityP() const
{
    return (velocityP);
}

/*! \brief Set P-wave velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP)
{
    velocityP = setVelocityP;
}

/*! \brief Get const reference to S-wave velocity
 * 
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getVelocityS()
{
    return (velocityS);
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getVelocityS() const
{
    return (velocityS);
}

/*! \brief Set S-wave velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS)
{
    velocityS = setVelocityS;
}

/*! \brief Get const reference to tauP
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getTauP()
{
    return (tauP);
}

/*! \brief Get const reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getTauP() const
{
    return (tauP);
}

/*! \brief Set tauP parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setTauP(scai::lama::Vector<ValueType> const &setTauP)
{
    tauP = setTauP;
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getTauS()
{
    return (tauS);
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getTauS() const
{
    return (tauS);
}

/*! \brief Set tauS parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setTauS(scai::lama::Vector<ValueType> const &setTauS)
{
    tauS = setTauS;
}

/*! \brief Get const reference to porosity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getPorosity()
{
    return (porosity);
}

/*! \brief Get const reference to porosity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getPorosity() const
{
    return (porosity);
}

/*! \brief Set porosity model parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setPorosity(scai::lama::Vector<ValueType> const &setPorosity)
{
    porosity = setPorosity;
}

/*! \brief Get const reference to saturation
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getSaturation()
{
    return (saturation);
}

/*! \brief Get const reference to saturation
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getSaturation() const
{
    return (saturation);
}

/*! \brief Set saturation model parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setSaturation(scai::lama::Vector<ValueType> const &setSaturation)
{
    saturation = setSaturation;
}

/*! \brief Get parameter normalizeGradient
 */
template <typename ValueType>
bool KITGPI::Gradient::Gradient<ValueType>::getNormalizeGradient() const
{
    return (normalizeGradient);
}

/*! \brief Set normalizeGradient parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setNormalizeGradient(bool const &gradNorm)
{
    normalizeGradient = gradNorm;
}

/*! \brief Overloading = Operation
 *
 \param rhs Gradient which is copied.
 */
template <typename ValueType>
KITGPI::Gradient::Gradient<ValueType> &KITGPI::Gradient::Gradient<ValueType>::operator=(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    assign(rhs);
    return *this;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Gradient<ValueType> &KITGPI::Gradient::Gradient<ValueType>::operator-=(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    minusAssign(rhs);
    return *this;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Gradient<ValueType> &KITGPI::Gradient::Gradient<ValueType>::operator+=(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    plusAssign(rhs);
    return *this;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Gradient<ValueType> &KITGPI::Gradient::Gradient<ValueType>::operator*=(ValueType const &rhs)
{
    timesAssign(rhs);
    return *this;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Gradient<ValueType> &KITGPI::Gradient::Gradient<ValueType>::operator*=(scai::lama::Vector<ValueType> const &rhs)
{
    timesAssign(rhs);
    return *this;
}

template class KITGPI::Gradient::Gradient<float>;
template class KITGPI::Gradient::Gradient<double>;
