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
    for (IndexType iComponent = 0; iComponent < KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (seismograms.getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            lama::DenseMatrix<ValueType> const &seismoA = seismograms.getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            init(seismoA.getRowDistributionPtr(), seismoA.getColDistributionPtr(), seismoA.getContextPtr());
        }
    }
}

/*! \brief Wrapper to support SeismogramHandler
 \param seismograms SeismogramHandler object
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::apply(KITGPI::Acquisition::SeismogramHandler<ValueType> &seismograms) const
{
    for (IndexType iComponent = 0; iComponent < KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
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

/*! \brief Initialize taper
 \param rowDist Row distribution
 \param colDist Column distribution
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::initModelTransform(dmemo::DistributionPtr dist, dmemo::DistributionPtr distEM, hmemo::ContextPtr ctx)
{
    modelTransformMatrixToSeismic.allocate(dist, distEM);
    modelTransformMatrixToSeismic.setContextPtr(ctx);
    modelTransformMatrixToEM.allocate(distEM, dist);
    modelTransformMatrixToEM.setContextPtr(ctx);
    modelParameterTransform.allocate(dist);
    modelParameterTransform.setContextPtr(ctx);
    modelParameterTransformEM.allocate(distEM);
    modelParameterTransformEM.setContextPtr(ctx);
}

/*! \brief Apply model transform to a parameter
 \param gradientParameter gradient parameter
 */
template <typename ValueType>
lama::Vector<ValueType> const &KITGPI::Taper::Taper2D<ValueType>::applyGradientTransformToEM(lama::Vector<ValueType> const &gradientParameter)
{
    modelParameterTransformEM = modelTransformMatrixToEM * gradientParameter;
        
    return modelParameterTransformEM;
}

/*! \brief Apply model transform to a parameter
 \param modelParameter model parameter
 */
template <typename ValueType>
lama::Vector<ValueType> const &KITGPI::Taper::Taper2D<ValueType>::applyModelTransformToEM(lama::Vector<ValueType> const &modelParameter, lama::Vector<ValueType> const &modelParameterEM)
{
    lama::DenseVector<ValueType> modelParameterResidualEM;
    
    modelParameterTransform = modelTransformMatrixToSeismic * modelParameterEM;
    modelParameterResidualEM = modelTransformMatrixToEM * modelParameterTransform;
    modelParameterResidualEM = modelParameterEM - modelParameterResidualEM;
    
    modelParameterTransformEM = modelTransformMatrixToEM * modelParameter;
    modelParameterTransformEM += modelParameterResidualEM;
        
    return modelParameterTransformEM;
}

/*! \brief Apply model transform to a parameter
 \param modelParameter model parameter
 */
template <typename ValueType>
lama::Vector<ValueType> const &KITGPI::Taper::Taper2D<ValueType>::applyGradientTransformToSeismic(lama::Vector<ValueType> const &gradientParameterEM)
{
    modelParameterTransform = modelTransformMatrixToSeismic * gradientParameterEM;
    
    return modelParameterTransform;
}

/*! \brief Apply model transform to a parameter
 \param modelParameter model parameter
 */
template <typename ValueType>
lama::Vector<ValueType> const &KITGPI::Taper::Taper2D<ValueType>::applyModelTransformToSeismic(lama::Vector<ValueType> const &modelParameter, lama::Vector<ValueType> const &modelParameterEM)
{
    lama::DenseVector<ValueType> modelParameterResidual;
    
    modelParameterTransformEM = modelTransformMatrixToEM * modelParameter;
    modelParameterResidual = modelTransformMatrixToSeismic * modelParameterTransformEM;
    modelParameterResidual = modelParameter - modelParameterResidual;
    
    modelParameterTransform = modelTransformMatrixToSeismic * modelParameterEM;
    modelParameterTransform += modelParameterResidual;
    
    return modelParameterTransform;
}

/*! \brief Read a taper from file
 * \param filename taper filename
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::read(std::string filename)
{
    dmemo::DistributionPtr distTraces(data.getRowDistributionPtr());
    dmemo::DistributionPtr distSamples(data.getColDistributionPtr());

    data.readFromFile(filename);

    data.redistribute(distTraces, distSamples);
}

/*! \brief calculate a matrix to transform Seismic model to EM model
 * \param modelCoordinates coordinates of the seismic model
 * \param config Configuration seismic
 * \param modelCoordinatesEM coordinates of the EM model
 * \param configEM Configuration EM
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::calcSeismictoEMMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesEM, KITGPI::Configuration::Configuration configEM)
{
    // this function is only valid for grid partitioning with useVariableGrid=0
    KITGPI::Acquisition::coordinate3D coordinate;
    KITGPI::Acquisition::coordinate3D coordinateEM;
    ValueType DH = modelCoordinates.getDH();
    ValueType DHEM = modelCoordinatesEM.getDH();
    ValueType x;
    ValueType y;
    ValueType z;
    ValueType xEM;
    ValueType yEM;
    ValueType zEM;
    ValueType x0 = modelCoordinates.getX0();
    ValueType y0 = modelCoordinates.getY0();
    ValueType z0 = modelCoordinates.getZ0();
    ValueType x0EM = modelCoordinatesEM.getX0();
    ValueType y0EM = modelCoordinatesEM.getY0();
    ValueType z0EM = modelCoordinatesEM.getZ0();
    IndexType supplementaryOrder;
    IndexType supplementaryWidth;        
    
    if (DH <= DHEM) {
        supplementaryOrder = 1;
    } else {        
        supplementaryOrder = std::round(DH/DHEM);
    }
    supplementaryWidth = std::floor((supplementaryOrder-1)/2);
    
    dmemo::DistributionPtr distEM(modelTransformMatrixToEM.getRowDistributionPtr());
    dmemo::DistributionPtr dist(modelTransformMatrixToEM.getColDistributionPtr());
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    IndexType rowIndex;
    lama::MatrixAssembly<ValueType> assembly;

    if (DH <= DHEM) {
        SparseFormat EMtoSeismicMatrix;
        EMtoSeismicMatrix.allocate(dist, distEM);
        distEM->getOwnedIndexes(ownedIndexes);
        for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
            coordinateEM = modelCoordinatesEM.index2coordinate(ownedIndex);
            xEM=coordinateEM.x*DHEM+x0EM;
            yEM=coordinateEM.y*DHEM+y0EM;
            zEM=coordinateEM.z*DHEM+z0EM;
            coordinate.x=round((xEM-x0)/DH)-supplementaryWidth;
            coordinate.y=round((yEM-y0)/DH)-supplementaryWidth;
            coordinate.z=round((zEM-z0)/DH)-supplementaryWidth;
            for (IndexType iXorder = 0; iXorder < supplementaryOrder; iXorder++) {            
                for (IndexType iYorder = 0; iYorder < supplementaryOrder; iYorder++) {            
                    for (IndexType iZorder = 0; iZorder < supplementaryOrder; iZorder++) {
                        if (coordinate.x+iXorder<modelCoordinates.getNX() && coordinate.y+iYorder<modelCoordinates.getNY() && coordinate.z+iZorder<modelCoordinates.getNZ() && coordinate.x+iXorder>=0 && coordinate.y+iYorder>=0 && coordinate.z+iZorder>=0) {
                            rowIndex=modelCoordinates.coordinate2index(coordinate.x+iXorder, coordinate.y+iYorder, coordinate.z+iZorder);
                            assembly.push(rowIndex, ownedIndex, 1);
                        }
                    }
                }
            }
        }
        EMtoSeismicMatrix.fillFromAssembly(assembly);
        modelTransformMatrixToEM = transpose(EMtoSeismicMatrix);
    } else {
        dist->getOwnedIndexes(ownedIndexes);
        for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
            coordinate = modelCoordinates.index2coordinate(ownedIndex);
            x=coordinate.x*DH+x0;
            y=coordinate.y*DH+y0;
            z=coordinate.z*DH+z0;
            coordinateEM.x=round((x-x0EM)/DHEM)-supplementaryWidth;
            coordinateEM.y=round((y-y0EM)/DHEM)-supplementaryWidth;
            coordinateEM.z=round((z-z0EM)/DHEM)-supplementaryWidth;
            for (IndexType iXorder = 0; iXorder < supplementaryOrder; iXorder++) {            
                for (IndexType iYorder = 0; iYorder < supplementaryOrder; iYorder++) {            
                    for (IndexType iZorder = 0; iZorder < supplementaryOrder; iZorder++) {
                        if (coordinateEM.x+iXorder<modelCoordinatesEM.getNX()  && coordinateEM.y+iYorder<modelCoordinatesEM.getNY()  && coordinateEM.z+iZorder<modelCoordinatesEM.getNZ()  && coordinateEM.x+iXorder>=0 && coordinateEM.y+iYorder>=0 && coordinateEM.z+iZorder>=0) {
                            rowIndex=modelCoordinatesEM.coordinate2index(coordinateEM.x+iXorder, coordinateEM.y+iYorder, coordinateEM.z+iZorder);
                            assembly.push(rowIndex, ownedIndex, 1);
                        }
                    }
                }
            }
        }
        modelTransformMatrixToEM.fillFromAssembly(assembly);
    }
//     modelTransformMatrixToEM.writeToFile("model/modelTransformMatrixToEM.mtx");
}

/*! \brief calculate a matrix to transform EM model to Seismic model
 * \param modelCoordinates coordinates of the seismic model
 * \param config Configuration seismic
 * \param modelCoordinatesEM coordinates of the EM model
 * \param configEM Configuration EM
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::calcEMtoSeismicMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesEM, KITGPI::Configuration::Configuration configEM)
{
    // this function is only valid for grid partitioning with useVariableGrid=0
    KITGPI::Acquisition::coordinate3D coordinate;
    KITGPI::Acquisition::coordinate3D coordinateEM;
    ValueType DH = modelCoordinates.getDH();
    ValueType DHEM = modelCoordinatesEM.getDH();
    ValueType x;
    ValueType y;
    ValueType z;
    ValueType xEM;
    ValueType yEM;
    ValueType zEM;
    ValueType x0 = modelCoordinates.getX0();
    ValueType y0 = modelCoordinates.getY0();
    ValueType z0 = modelCoordinates.getZ0();
    ValueType x0EM = modelCoordinatesEM.getX0();
    ValueType y0EM = modelCoordinatesEM.getY0();
    ValueType z0EM = modelCoordinatesEM.getZ0();
    IndexType supplementaryOrder;
    IndexType supplementaryWidth;        
    
    if (DHEM <= DH) {
        supplementaryOrder = 1;
    } else {        
        supplementaryOrder = std::round(DHEM/DH);
    }
    supplementaryWidth = std::floor((supplementaryOrder-1)/2);
    
    dmemo::DistributionPtr dist(modelTransformMatrixToSeismic.getRowDistributionPtr());
    dmemo::DistributionPtr distEM(modelTransformMatrixToSeismic.getColDistributionPtr());
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    IndexType rowIndex;
    lama::MatrixAssembly<ValueType> assembly;

    if (DH <= DHEM) {
        distEM->getOwnedIndexes(ownedIndexes);
        for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
            coordinateEM = modelCoordinatesEM.index2coordinate(ownedIndex);
            xEM=coordinateEM.x*DHEM+x0EM;
            yEM=coordinateEM.y*DHEM+y0EM;
            zEM=coordinateEM.z*DHEM+z0EM;
            coordinate.x=round((xEM-x0)/DH)-supplementaryWidth;
            coordinate.y=round((yEM-y0)/DH)-supplementaryWidth;
            coordinate.z=round((zEM-z0)/DH)-supplementaryWidth;
            for (IndexType iXorder = 0; iXorder < supplementaryOrder; iXorder++) {            
                for (IndexType iYorder = 0; iYorder < supplementaryOrder; iYorder++) {            
                    for (IndexType iZorder = 0; iZorder < supplementaryOrder; iZorder++) {
                        if (coordinate.x+iXorder<modelCoordinates.getNX() && coordinate.y+iYorder<modelCoordinates.getNY() && coordinate.z+iZorder<modelCoordinates.getNZ() && coordinate.x+iXorder>=0 && coordinate.y+iYorder>=0 && coordinate.z+iZorder>=0) {
                            rowIndex=modelCoordinates.coordinate2index(coordinate.x+iXorder, coordinate.y+iYorder, coordinate.z+iZorder);
                            assembly.push(rowIndex, ownedIndex, 1);
                        }
                    }
                }
            }
        }
        modelTransformMatrixToSeismic.fillFromAssembly(assembly);
    } else {
        SparseFormat SeismictoEMMatrix;
        SeismictoEMMatrix.allocate(distEM, dist);
        dist->getOwnedIndexes(ownedIndexes);
        for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
            coordinate = modelCoordinates.index2coordinate(ownedIndex);
            x=coordinate.x*DH+x0;
            y=coordinate.y*DH+y0;
            z=coordinate.z*DH+z0;
            coordinateEM.x=round((x-x0EM)/DHEM)-supplementaryWidth;
            coordinateEM.y=round((y-y0EM)/DHEM)-supplementaryWidth;
            coordinateEM.z=round((z-z0EM)/DHEM)-supplementaryWidth;
            for (IndexType iXorder = 0; iXorder < supplementaryOrder; iXorder++) {            
                for (IndexType iYorder = 0; iYorder < supplementaryOrder; iYorder++) {            
                    for (IndexType iZorder = 0; iZorder < supplementaryOrder; iZorder++) {
                        if (coordinateEM.x+iXorder<modelCoordinatesEM.getNX()  && coordinateEM.y+iYorder<modelCoordinatesEM.getNY()  && coordinateEM.z+iZorder<modelCoordinatesEM.getNZ()  && coordinateEM.x+iXorder>=0 && coordinateEM.y+iYorder>=0 && coordinateEM.z+iZorder>=0) {
                            rowIndex=modelCoordinatesEM.coordinate2index(coordinateEM.x+iXorder, coordinateEM.y+iYorder, coordinateEM.z+iZorder);
                            assembly.push(rowIndex, ownedIndex, 1);
                        }
                    }
                }
            }
        }
        SeismictoEMMatrix.fillFromAssembly(assembly);
        modelTransformMatrixToSeismic = transpose(SeismictoEMMatrix);
    }
//     modelTransformMatrixToSeismic.writeToFile("model/modelTransformMatrixToSeismic.mtx");
}

/*! \brief exchange porosity and saturation from EM to seismic or from seismic to EM
 * \param model Seismic model
 * \param modelEM EM model
 * \param config EM config
 */
template <typename ValueType> void KITGPI::Taper::Taper2D<ValueType>::exchangePetrophysics(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Modelparameter::Modelparameter<ValueType> &modelEM, KITGPI::Configuration::Configuration config)
{
    lama::DenseVector<ValueType> porositytemp;
    lama::DenseVector<ValueType> saturationtemp;
    
    std::string equationType = config.get<std::string>("equationType");
    bool isSeismic = Common::checkEquationType<ValueType>(equationType); 
    IndexType exchangeStrategy = config.get<IndexType>("exchangeStrategy");  
    if (!isSeismic) {        
        if (exchangeStrategy == 2) {
            // case 2: exchange all of the petrophysical parameters 
            porositytemp = this->applyModelTransformToEM(model.getPorosity(), modelEM.getPorosity()); 
            saturationtemp = this->applyModelTransformToEM(model.getSaturation(), modelEM.getSaturation()); 
                    
            modelEM.setPorosity(porositytemp);
            modelEM.setSaturation(saturationtemp);           
        } else if (exchangeStrategy != 0) {
            // case 1,3,4,5,6: exchange porosity to EM model  
            porositytemp = this->applyModelTransformToEM(model.getPorosity(), modelEM.getPorosity());
                    
            modelEM.setPorosity(porositytemp);
        }  
        // case 0,1,2,3,4: self-constraint of the petrophysical relationship   
        modelEM.calcWaveModulusFromPetrophysics();   
        if (config.get<bool>("useModelThresholds"))
            modelEM.applyThresholds(config); 
    } else {
        if (exchangeStrategy == 2) {
            // case 2: exchange all of the petrophysical parameters 
            porositytemp = this->applyModelTransformToSeismic(model.getPorosity(), modelEM.getPorosity()); 
            saturationtemp = this->applyModelTransformToSeismic(model.getSaturation(), modelEM.getSaturation()); 
                    
            model.setPorosity(porositytemp);
            model.setSaturation(saturationtemp);           
        } else if (exchangeStrategy != 0) {     
            // case 1,3,4,5,6: exchange saturation to seismic model             
            saturationtemp = this->applyModelTransformToSeismic(model.getSaturation(), modelEM.getSaturation()); 
                    
            model.setSaturation(saturationtemp); 
        }
        // case 0,1,2,3,4: self-constraint of the petrophysical relationship                     
        model.calcWaveModulusFromPetrophysics(); 
        if (config.get<bool>("useModelThresholds"))
            model.applyThresholds(config); 
    }
}

/*! \brief exchange mode parameters from one seismic/EM wave to another seismic/EM wave
 * \param model1 Seismic/EM model
 * \param model2 Seismic/EM model
 * \param config1 Seismic/EM config
 * \param config2 Seismic/EM config
 */
template <typename ValueType> void KITGPI::Taper::Taper2D<ValueType>::exchangeModelparameters(KITGPI::Modelparameter::Modelparameter<ValueType> &model1, KITGPI::Modelparameter::Modelparameter<ValueType> &model2, KITGPI::Configuration::Configuration config1, KITGPI::Configuration::Configuration config2)
{
    std::string equationType1 = config1.get<std::string>("equationType");
    std::transform(equationType1.begin(), equationType1.end(), equationType1.begin(), ::tolower);  
    bool isSeismic1 = Common::checkEquationType<ValueType>(equationType1); 
    IndexType exchangeStrategy1 = config1.get<IndexType>("exchangeStrategy");  
    std::string equationType2 = config2.get<std::string>("equationType");
    std::transform(equationType2.begin(), equationType2.end(), equationType2.begin(), ::tolower);  
    bool isSeismic2 = Common::checkEquationType<ValueType>(equationType2); 
    IndexType exchangeStrategy2 = config2.get<IndexType>("exchangeStrategy"); 
    SCAI_ASSERT_ERROR(isSeismic1 == isSeismic2, "isSeismic1 != isSeismic2");
    SCAI_ASSERT_ERROR(exchangeStrategy1 == exchangeStrategy2, "exchangeStrategy1 != exchangeStrategy2");
    SCAI_ASSERT_ERROR((exchangeStrategy1 > 1), "exchangeStrategy mush be 0, 2, 4, 6 in exchangeModelparameters()");
    
    model2.setPorosity(model1.getPorosity());
    model2.setSaturation(model1.getSaturation());
    model2.setReflectivity(model1.getReflectivity());
    if (isSeismic1 && exchangeStrategy1 > 1) {        
        model2.setDensity(model1.getDensity());
        if ((equationType1.compare("sh") == 0 || equationType1.compare("viscosh") == 0 || equationType1.compare("elastic") == 0 || equationType1.compare("viscoelastic") == 0) && (equationType2.compare("sh") == 0 || equationType2.compare("viscosh") == 0 || equationType2.compare("elastic") == 0 || equationType2.compare("viscoelastic") == 0)) {
            model2.setVelocityS(model1.getVelocityS());
        }
        if ((equationType1.compare("acoustic") == 0 || equationType1.compare("elastic") == 0 || equationType1.compare("viscoelastic") == 0) && (equationType2.compare("acoustic") == 0 || equationType2.compare("elastic") == 0 || equationType2.compare("viscoelastic") == 0)) {
            model2.setVelocityP(model1.getVelocityP());
        }
        if ((equationType1.compare("viscoelastic") == 0 || equationType1.compare("viscosh") == 0) && (equationType2.compare("viscoelastic") == 0 || equationType2.compare("viscosh") == 0)){
            model2.setTauS(model1.getTauS());
        }
        if (equationType1.compare("viscoelastic") == 0 && equationType2.compare("viscoelastic") == 0){
            model2.setTauP(model1.getTauP());
        }
    } else if (!isSeismic1 && exchangeStrategy1 > 1) {
        model2.setMagneticPermeability(model1.getMagneticPermeability());
        model2.setElectricConductivity(model1.getElectricConductivity());
        model2.setDielectricPermittivity(model1.getDielectricPermittivity());
        if ((equationType1.compare("viscotmem") == 0 || equationType1.compare("viscoemem") == 0) && (equationType2.compare("viscotmem") == 0 || equationType2.compare("viscoemem") == 0)){
            model2.setTauElectricConductivity(model1.getTauElectricConductivity());
            model2.setTauDielectricPermittivity(model1.getTauDielectricPermittivity());
        }
    }     
}

/*! \brief Initialize wavefield transform matrix
 \param config Configuration
 \param distInversion distribution inversion
 \param dist distribution
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::initAverageMatrix(KITGPI::Configuration::Configuration config, dmemo::DistributionPtr distInversion, dmemo::DistributionPtr dist, hmemo::ContextPtr ctx)
{
    averageMatrix = lama::zero<SparseFormat>(distInversion, dist);
    averageMatrix.setContextPtr(ctx);
    recoverMatrix = lama::zero<SparseFormat>(dist, distInversion);
    recoverMatrix.setContextPtr(ctx);
}

/*! \brief calculate a matrix for averaging inversion
 * \param modelCoordinates coordinates of the original model
 * \param modelCoordinatesInversion coordinates of the averaged model
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::calcAverageMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesInversion)
{
    // this function is only valid for grid partitioning with useVariableGrid=0
    ValueType averageValue;
    IndexType NX = modelCoordinates.getNX();
    IndexType NY = modelCoordinates.getNY();
    IndexType NZ = modelCoordinates.getNZ();
    ValueType DHInversion = ceil(NX/modelCoordinatesInversion.getNX());
    KITGPI::Acquisition::coordinate3D coordinate;
    KITGPI::Acquisition::coordinate3D coordinateInversion;
    
    if (NZ>1) {
        averageValue=1/(pow(DHInversion,3));
    } else {
        averageValue=1/(pow(DHInversion,2));
    }
    
    dmemo::DistributionPtr distInversion(averageMatrix.getRowDistributionPtr());
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

    averageMatrix.fillFromAssembly(assemblyAverage);
    recoverMatrix.fillFromAssembly(assemblyRecover);
//     averageMatrix.writeToFile("model/averageMatrix_" + std::to_string(NX) + "_" + std::to_string(NY) + "_" + std::to_string(NZ) + "_" + std::to_string(averageValue) + ".mtx");
//     recoverMatrix.writeToFile("model/recoverMatrix_" + std::to_string(NX) + "_" + std::to_string(NY) + "_" + std::to_string(NZ) + ".mtx");
}

/*! \brief Get averageMatrix
 */
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::Taper::Taper2D<ValueType>::getAverageMatrix()
{
    return averageMatrix;
}

/*! \brief Get recoverMatrix
 */
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::Taper::Taper2D<ValueType>::getRecoverMatrix()
{
    return recoverMatrix;
}

template class KITGPI::Taper::Taper2D<double>;
template class KITGPI::Taper::Taper2D<float>;
