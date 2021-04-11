#include "Taper2D.hpp"

using namespace scai;

/*! \brief Initialize taper
 \param rowDist Row distribution
 \param colDist Column distribution
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::init(dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist, hmemo::ContextPtr ctx)
{
    data.allocate(rowDist, colDist);
    data.setContextPtr(ctx);
}

/*! \brief Initialize taper
 \param seismograms seismograms hander
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::init(KITGPI::Acquisition::SeismogramHandler<ValueType> const seismograms)
{
    for (scai::IndexType iComponent = 0; iComponent < 4; iComponent++) {
        if (seismograms.getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            lama::DenseMatrix<ValueType> const &seismoA = seismograms.getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            init(seismoA.getRowDistributionPtr(), seismoA.getColDistributionPtr(), seismoA.getContextPtr());
        }
    }
}

/*! \brief Initialize taper
 \param rowDist Row distribution
 \param colDist Column distribution
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::initWavefieldTransform(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr averageDist, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, bool isSeismic)
{
    wavefieldAverageMatrix = scai::lama::zero<SparseFormat>(averageDist, dist);
    wavefieldAverageMatrix.setContextPtr(ctx);
    wavefieldRecoverMatrix = scai::lama::zero<SparseFormat>(dist, averageDist);
    wavefieldRecoverMatrix.setContextPtr(ctx);
    
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);   
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);  
    if (isSeismic) {
        wavefieldAverage = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
        wavefieldRecover = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
        wavefieldAverage->init(ctx, averageDist);
        wavefieldRecover->init(ctx, dist);   
    }
}

/*! \brief Wrapper to support SeismogramHandler
 \param seismograms SeismogramHandler object
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::apply(KITGPI::Acquisition::SeismogramHandler<ValueType> &seismograms) const
{
    for (scai::IndexType iComponent = 0; iComponent < 4; iComponent++) {
        if (seismograms.getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            Acquisition::Seismogram<ValueType> &thisSeismogram = seismograms.getSeismogram(Acquisition::SeismogramType(iComponent));
            apply(thisSeismogram);
        }
    }
}


/*! \brief Apply taper to a single seismogram
 \param seismogram Seismogram
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const
{
    lama::DenseMatrix<ValueType> &seismogramData = seismogram.getData();
    seismogramData.binaryOp(seismogramData, common::BinaryOp::MULT, data);
}

/*! \brief Apply taper to a DenseMatrix
 \param mat Seismogram
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::apply(lama::DenseMatrix<ValueType> &mat) const
{
    mat.binaryOp(mat, common::BinaryOp::MULT, data);
}

/*! \brief Apply model transform to a wavefield
 \param modelParameter model parameter
 */
template <typename ValueType>
KITGPI::Wavefields::Wavefields<ValueType> &KITGPI::Taper::Taper2D<ValueType>::applyWavefieldAverage(typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr &wavefieldPtr)
{
    wavefieldAverage->applyWavefieldTransform(wavefieldAverageMatrix, *wavefieldPtr);
    return *wavefieldAverage;
}

/*! \brief Apply model recover to a wavefield
 \param modelParameter model parameter
 */
template <typename ValueType>
KITGPI::Wavefields::Wavefields<ValueType> &KITGPI::Taper::Taper2D<ValueType>::applyWavefieldRecover(typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr &wavefieldPtr)
{
    wavefieldRecover->applyWavefieldTransform(wavefieldRecoverMatrix, *wavefieldPtr);
    return *wavefieldRecover;
}

/*! \brief Read a taper from file
 * \param filename taper filename
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::read(std::string filename)
{
    scai::dmemo::DistributionPtr distTraces(data.getRowDistributionPtr());
    scai::dmemo::DistributionPtr distSamples(data.getColDistributionPtr());

    data.readFromFile(filename);

    data.redistribute(distTraces, distSamples);
}

/*! \brief calculate an average matrix for inversion
 * \param modelCoordinates coordinates of the original model
 * \param modelCoordinatesInversion coordinates of the averaged model
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::calcInversionAverageMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesInversion)
{
    // this function is only valid for grid partitioning with useVariableGrid=0
    ValueType averageValue;
    scai::IndexType NX = modelCoordinates.getNX();
    scai::IndexType NY = modelCoordinates.getNY();
    scai::IndexType NZ = modelCoordinates.getNZ();
    ValueType DHInversion = ceil(NX/modelCoordinatesInversion.getNX());
    KITGPI::Acquisition::coordinate3D coordinate;
    KITGPI::Acquisition::coordinate3D coordinateInversion;
    
    if (NZ>1) {
        averageValue=1/(pow(DHInversion,3));
    } else {
        averageValue=1/(pow(DHInversion,2));
    }
    
    scai::dmemo::DistributionPtr distInversion(wavefieldAverageMatrix.getRowDistributionPtr());
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    distInversion->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assemblyAverage;
    lama::MatrixAssembly<ValueType> assemblyRecover;
    IndexType columnIndex;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
        coordinateInversion = modelCoordinatesInversion.index2coordinate(ownedIndex);
        coordinate.x = coordinateInversion.x*DHInversion;
        coordinate.y = coordinateInversion.y*DHInversion;
        coordinate.z = coordinateInversion.z*DHInversion;

        for (IndexType iXorder = 0; iXorder < DHInversion; iXorder++) {            
            for (IndexType iYorder = 0; iYorder < DHInversion; iYorder++) {            
                for (IndexType iZorder = 0; iZorder < DHInversion; iZorder++) {
                    if ((coordinate.x+iXorder < NX) && (coordinate.y+iYorder < NY) && (coordinate.z+iZorder < NZ)) {
                        columnIndex=modelCoordinates.coordinate2index(coordinate.x+iXorder, coordinate.y+iYorder, coordinate.z+iZorder);
                        assemblyAverage.push(ownedIndex, columnIndex, averageValue);
                        assemblyRecover.push(columnIndex, ownedIndex, 1);
                    }
                }
            }
        }
    }

    wavefieldAverageMatrix.fillFromAssembly(assemblyAverage);
    wavefieldRecoverMatrix.fillFromAssembly(assemblyRecover);
    wavefieldAverageMatrix.writeToFile("model/wavefieldAverageMatrix_" + std::to_string(NX) + "_" + std::to_string(NY) + "_" + std::to_string(NZ) + "_" + std::to_string(averageValue) + ".mtx");
    wavefieldRecoverMatrix.writeToFile("model/wavefieldRecoverMatrix_" + std::to_string(NX) + "_" + std::to_string(NY) + "_" + std::to_string(NZ) + ".mtx");
}

template class KITGPI::Taper::Taper2D<double>;
template class KITGPI::Taper::Taper2D<float>;
