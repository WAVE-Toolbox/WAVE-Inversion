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
}

/*! \brief If stream configuration is used, calculate a weighting vector to balance gradient per shot
 \param modelPerShot model per shot
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::calcWeightingVector(scai::lama::Vector<ValueType> const &gradientPerShotVector, scai::IndexType NY, scai::IndexType shotInd)
{
    auto dist = gradientPerShotVector.getDistributionPtr();
    scai::hmemo::ContextPtr ctx = gradientPerShotVector.getContextPtr();

    scai::lama::DenseVector<ValueType> weightingVector;
    scai::lama::DenseVector<ValueType> gradientPerShotVectorTemp;
    weightingVector.allocate(dist);
    weightingVector.setContextPtr(ctx);
    gradientPerShotVectorTemp.unaryOp(gradientPerShotVector, common::UnaryOp::ABS); 
    if (weightGradient == 0) {
        weightingVector = 1.0;
    } else if (weightGradient == 1) {
        weightingVector = 0.0;
        IndexType NX = gradientPerShotVectorTemp.size() / NY;
        for (IndexType ix = 0; ix < NX; ix++) {
            for (IndexType iy = 0; iy < NY; iy++) {
                weightingVector[ix] = weightingVector[ix] + gradientPerShotVectorTemp[iy*NX+ix];
            }
            for (IndexType iy = 1; iy < NY; iy++) {
                weightingVector[iy*NX+ix] = weightingVector[ix];
            }
        }
        if (weightingVector.maxNorm() != 0)
            weightingVector *= 1.0 / weightingVector.maxNorm();
    }
//     IO::writeVector(weightingVector, "gradients/weightingVector.shot_" + std::to_string(shotInd+1), 1);
    
    return weightingVector;
}

/*! \brief get an inner workflow */
template <typename ValueType>
std::vector<bool> KITGPI::Gradient::Gradient<ValueType>::getInvertForParameters()
{
    return workflowInner.getInvertForParameters();
}

/*! \brief get an inner workflow */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::printInvertForParameters(scai::dmemo::CommunicatorPtr comm)
{
    workflowInner.printInvertForParameters(comm);
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

/*! \brief Get const reference to reflectivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getReflectivity()
{
    return (reflectivity);
}

/*! \brief Get const reference to reflectivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getReflectivity() const
{
    return (reflectivity);
}

/*! \brief Set reflectivity
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setReflectivity(scai::lama::Vector<ValueType> const &setReflectivity)
{
    reflectivity = setReflectivity;
}

/*! \brief Get parameter normalizeGradient
 */
template <typename ValueType>
bool KITGPI::Gradient::Gradient<ValueType>::getNormalizeGradient() const
{
    return (normalizeGradient);
}

/*! \brief prepare for inversion
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::prepareForInversion(KITGPI::Configuration::Configuration config)
{
    normalizeGradient = config.get<bool>("normalizeGradient");
    weightGradient = config.getAndCatch("weightGradient", 0);
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
