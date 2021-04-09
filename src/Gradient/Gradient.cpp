#include "Gradient.hpp"
#include <IO/IO.hpp>

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
    IO::writeVector(vector,filename,fileFormat);
};

/*! \brief Read a parameter from file
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::readParameterisation(scai::lama::Vector<ValueType> &vector, std::string filename, IndexType fileFormat)
{    
    IO::readVector(vector,filename,fileFormat);
};

/*! \brief Allocate a single parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::allocateParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

/*! \brief Get const reference to density  parameter
 * 
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Gradient<ValueType>::getDensity()
{

    return (density);
}

/*! \brief Get const reference to density  parameter
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

/*! \brief Set tauP velocity  parameter
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

/*! \brief Get parameter normalizeGradient
 */
template <typename ValueType>
bool KITGPI::Gradient::Gradient<ValueType>::getNormalizeGradient() const
{
    return (normalizeGradient);
}

/*! \brief Set tauS velocity  parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setTauS(scai::lama::Vector<ValueType> const &setTauS)
{
    tauS = setTauS;
}

/*! \brief Set normalizeGradient parameter
 */
template <typename ValueType>
void KITGPI::Gradient::Gradient<ValueType>::setNormalizeGradient(bool const &gradNorm)
{
    normalizeGradient = gradNorm;
}

/*! \brief Get matrix that multiplies with model matrices to get a subset
 \param dist Distribution of the subset
 \param distBig Distribution of the big model
 \param modelCoordinates coordinate class object of the subset
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinates coordinate where to cut the subset
 */
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> KITGPI::Gradient::Gradient<ValueType>::getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D &cutCoordinates)
{
    SparseFormat shrinkMatrix; //!< Shrink Multiplication matrix
    shrinkMatrix.allocate(dist,distBig);
      
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);
    lama::MatrixAssembly<ValueType> assembly;
    IndexType indexBig = 0;
 
    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) { // loop over all indices
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex); // get submodel coordinate from singleIndex
        coordinate.x = coordinate.x + cutCoordinates.x; // offset depends on shot number
        indexBig = modelCoordinatesBig.coordinate2index(coordinate);
        assembly.push(ownedIndex, indexBig, 1.0);
    }
    shrinkMatrix.fillFromAssembly(assembly);
    return shrinkMatrix;
}

/*! \brief Get erase-matrix that erases the old values in the big model
 \param dist Distribution of the subset
 \param distBig Distribution of the big model
 \param modelCoordinates coordinate class object of the subset
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinates coordinate where to cut the subset
 */
template <typename ValueType>
scai::lama::SparseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::getEraseVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D &cutCoordinates, scai::IndexType NX, scai::IndexType NYBig, scai::IndexType boundaryWidth)
{
    scai::lama::SparseVector<ValueType> eraseVector(distBig,1); //!< Shrink Multiplication matrix
      
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);
    lama::VectorAssembly<ValueType> assembly;
    IndexType indexBig = 0;
 
    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) { // loop over all indices
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex); // get submodel coordinate from singleIndex
        coordinate.x = coordinate.x + cutCoordinates.x; // offset depends on shot number
        indexBig = modelCoordinatesBig.coordinate2index(coordinate);
        assembly.push(indexBig, 0);
    }
    eraseVector.fillFromAssembly(assembly);
    
    // damp the boundary boarders
    for (IndexType y = 0; y < NYBig; y++) {
        for (IndexType i = 0; i < boundaryWidth; i++) {
            double tmp = (double)1.0-(i+1) / (double)boundaryWidth;
            eraseVector[modelCoordinatesBig.coordinate2index(cutCoordinates.x+i, y, 0)] = tmp;
            eraseVector[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-1-i, y, 0)] = tmp;
        }
    }
    return eraseVector;
}

/*! \brief Smoothes parameter within a certain range via gauss
 \param modelCoordinatesBig coordinate class object of the big model
 \param parameter input parameter that is to be smoothed
 \param subsetSize size of the subset in x-direction
 \param cutCoordinates coordinate where to cut the subset
 \param smoothRange grid points to the left/right which are to be smoothed
 \param NX NX in model
 \param NXBig NX in big model
 \param NYBig NY in big model
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Gradient::Gradient<ValueType>::smoothParameter(Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, scai::lama::DenseVector<ValueType> parameter, Acquisition::coordinate3D &cutCoordinates, scai::IndexType smoothRange, scai::IndexType NX, scai::IndexType NXBig, scai::IndexType NYBig)
{
    scai::lama::DenseVector<ValueType> savedPar = parameter;
 
    for (IndexType y = 0; y < NYBig; y++) {
        for (IndexType i = 0; i<smoothRange*2+1; i++) {
            if (cutCoordinates.x < smoothRange+3) {
                parameter[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i, y, 0)] = savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i-3, y, 0)]*0.0055 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i-2, y, 0)]*0.061 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i-1, y, 0)]*0.242 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i, y, 0)]*0.383 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i+1, y, 0)]*0.242 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i+2, y, 0)]*0.061 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i+3, y, 0)]*0.0055;
            }
            else if (cutCoordinates.x+NX+smoothRange+3 > NXBig){
                parameter[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i, y, 0)] = savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i-3, y, 0)]*0.0055 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i-2, y, 0)]*0.061 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i-1, y, 0)]*0.242 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i, y, 0)]*0.383 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i+1, y, 0)]*0.242 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i+2, y, 0)]*0.061 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i+3, y, 0)]*0.0055;
            }
            else {
                parameter[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i, y, 0)] = savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i-3, y, 0)]*0.0055 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i-2, y, 0)]*0.061 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i-1, y, 0)]*0.242 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i, y, 0)]*0.383 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i+1, y, 0)]*0.242 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i+2, y, 0)]*0.061 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x-smoothRange+i+3, y, 0)]*0.0055;
            
                parameter[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i, y, 0)] = savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i-3, y, 0)]*0.0055 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i-2, y, 0)]*0.061 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i-1, y, 0)]*0.242 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i, y, 0)]*0.383 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i+1, y, 0)]*0.242 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i+2, y, 0)]*0.061 + savedPar[modelCoordinatesBig.coordinate2index(cutCoordinates.x+NX-smoothRange+i+3, y, 0)]*0.0055;
            }
        }
    }
    return parameter;
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
