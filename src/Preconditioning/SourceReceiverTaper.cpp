#include "SourceReceiverTaper.hpp"
using namespace scai;

/*! \brief Get taper
 *
 */
template <typename ValueType>
scai::lama::SparseVector<ValueType> KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::getTaper()
{
    return taper;
}

/*! \brief Get taper
 *
 */
template <typename ValueType>
std::vector<scai::lama::SparseVector<ValueType>> KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::getTaperEncode()
{
    return taperEncode;
}

/*! \brief stack taper 
 *
 \param gradient gradient
 */
template <typename ValueType>
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::stack(scai::lama::SparseVector<ValueType> taperSingle)
{   
    IndexType numshotsSuperShot = taperEncode.size();
    for (scai::IndexType shotInd = 0; shotInd < numshotsSuperShot; shotInd++) {
        taperEncode[shotInd] *= taperSingle;
    }
}

/*! \brief apply taper on gradient
 *
 \param gradient gradient
 */
template <typename ValueType>
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::apply(KITGPI::Gradient::Gradient<ValueType> &gradient)
{   
    gradient *= taper;
}

/*! \brief Init taper
 *
 \param dist_wavefield Distribution of the wavefields
 \param ctx Context
 \param Acquisition sources or receivers
 \param config Configuration class, which is used to derive all required parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param radius radius of the taper. 0=not apply taper.
 */
template <typename ValueType>
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Acquisition::AcquisitionGeometry<ValueType> const &Acquisition, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType radius)
{    
    /* Get local "global" indices */
    hmemo::HArray<IndexType> ownedIndexes;
    dist->getOwnedIndexes(ownedIndexes); // get global indexes for local part

    //gets coordinate index for first source to be changed
    IndexType SourceCoordinate;
    
    IndexType taperType = config.get<IndexType>("sourceReceiverTaperType");

    taper.setSameValue(dist, 1.0); // sparse vector taper gets zero element 1.0

    lama::VectorAssembly<ValueType> assembly;

    // distributed acquisition coordinates needed for each processor
    auto acquisition1DCoordinates = Acquisition.get1DCoordinates();
    acquisition1DCoordinates.replicate();
    auto rAcquisition1DCoordinates = scai::hmemo::hostReadAccess(acquisition1DCoordinates.getLocalValues());

    if(radius>0 && taperType>0){
        //loop over all gridpoint -> this is costly maybe there is a better alternative because the size of the taper is known
        for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)){
            ValueType valueTemp=1;
            ValueType value=1;
            
            KITGPI::Acquisition::coordinate3D coord;
            // get 3D coordinate of current ownedIndex
            coord = modelCoordinates.index2coordinate(ownedIndex);
                
            //loop over all shots or receivers
            for (IndexType srcRecNum = 0; srcRecNum < acquisition1DCoordinates.size(); srcRecNum++) {
                // get 1D coordinate
                SourceCoordinate = rAcquisition1DCoordinates[srcRecNum];

                // get 3D coordinate of source or receiver position
                KITGPI::Acquisition::coordinate3D centerCoord = modelCoordinates.index2coordinate(SourceCoordinate);
                
                ValueType distance;

                // distance between src/rec position and current local gridpoint
                distance = sqrt((centerCoord.x - coord.x) * (centerCoord.x - coord.x) + (centerCoord.y - coord.y) * (centerCoord.y - coord.y) + (centerCoord.z - coord.z) * (centerCoord.z - coord.z));
            
                // taper area:
                if (distance <= radius) {
                    if (distance >= 1) {
                        //-1 is needed because sparseVector=Densevector sets zero element to 0 for (grad * taper) the zero element should be 1 -> sparsevector+1  
                        switch (taperType) {
                            case 1:
                                //normalized logarithmic function : log(i)/log(max) for log(i)>1 
                                valueTemp=std::log(distance) / std::log(radius);
                                
                                break;
                            case 2:
                                //cos^2 function : cos(i*pi/(2*max))^2
                                valueTemp=common::Math::pow<ValueType>(common::Math::sin<ValueType>(distance * M_PI / (2.0 * radius)), 2.0);
                        }
                        
                        if (valueTemp < value){
                            value=valueTemp;
                        }
                    } else {
                        value=0;
                    }                    
                }
            }
            assembly.push(ownedIndex, value);
        }
    }

    taper.setContextPtr(ctx);
    taper.fillFromAssembly(assembly);
}

/*! \brief Init taper
 *
 \param dist_wavefield Distribution of the wavefields
 \param ctx Context
 \param sources sources
 \param receivers receivers
 \param config Configuration class, which is used to derive all required parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param radius radius of the taper. 0=not apply taper.
 */
template <typename ValueType>
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Acquisition::AcquisitionGeometry<ValueType> const &sources, KITGPI::Acquisition::AcquisitionGeometry<ValueType> const &receivers, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    /* Get local "global" indices */
    hmemo::HArray<IndexType> ownedIndexes;
    dist->getOwnedIndexes(ownedIndexes); // get global indexes for local part
    
    IndexType taperType = 0;
    
    taper.setSameValue(dist, 1.0); // sparse vector taper gets zero element 1.0

    lama::VectorAssembly<ValueType> assembly;

    // distributed acquisition coordinates needed for each processor
    auto sources1DCoordinates = sources.get1DCoordinates();
    sources1DCoordinates.replicate();
    auto rSources1DCoordinates = scai::hmemo::hostReadAccess(sources1DCoordinates.getLocalValues());
    auto receivers1DCoordinates = receivers.get1DCoordinates();
    receivers1DCoordinates.replicate();
    auto rReceivers1DCoordinates = scai::hmemo::hostReadAccess(receivers1DCoordinates.getLocalValues());

    KITGPI::Acquisition::coordinate3D sourceCoordFirst = modelCoordinates.index2coordinate(rSources1DCoordinates[0]);
    KITGPI::Acquisition::coordinate3D receiverCoordLast = modelCoordinates.index2coordinate(rReceivers1DCoordinates[rReceivers1DCoordinates.size()-1]);
    
    if (rReceivers1DCoordinates.size() > 4) {
        KITGPI::Acquisition::coordinate3D receiverCoordFirst = modelCoordinates.index2coordinate(rReceivers1DCoordinates[rReceivers1DCoordinates.size()-5]);
        if(abs(receiverCoordLast.x - receiverCoordFirst.x) < abs(receiverCoordLast.y - receiverCoordFirst.y)){ // borehole
            taperType = 1; 
        }
    }
    if(taperType > 0){
        //loop over all gridpoint -> this is costly maybe there is a better alternative because the size of the taper is known
        for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)){
            ValueType value=1;
            
            KITGPI::Acquisition::coordinate3D coord;
            // get 3D coordinate of current ownedIndex
            coord = modelCoordinates.index2coordinate(ownedIndex);
            
            if ((coord.x - sourceCoordFirst.x) * (coord.x - receiverCoordLast.x) >= 0) {
                value=0;
            }
            assembly.push(ownedIndex, value); 
        }
    }

    taper.setContextPtr(ctx);
    taper.fillFromAssembly(assembly);
    
    IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    if (useSourceEncode == 0) {
        taperEncode.clear();
        taperEncode.push_back(taper);
    }
}

/*! \brief Init taper
 *
 \param dist_wavefield Distribution of the wavefields
 \param ctx Context
 \param sources sources
 \param receivers receivers
 \param config Configuration class, which is used to derive all required parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param radius radius of the taper. 0=not apply taper.
 */
template <typename ValueType>
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::init(scai::dmemo::CommunicatorPtr commShot, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, scai::IndexType shotNumberEncode)
{
    IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    if (useSourceEncode != 0) {
        double start_t_shot, end_t_shot; /* For timing */
        start_t_shot = common::Walltime::get();
                
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings;  
        Acquisition::Sources<ValueType> sources;
        ValueType shotIncr = 0; 
        sources.getAcquisitionSettings(config, shotIncr); // to get numshots
        bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
        if (useStreamConfig) {
            Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
            Acquisition::Coordinates<ValueType> modelCoordinatesBig(configBig);
            std::vector<Acquisition::coordinate3D> cutCoordinates;
            std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;
            sourceSettingsBig = sources.getSourceSettings(); 
            Acquisition::getCutCoord(config, cutCoordinates, sourceSettingsBig, modelCoordinates, modelCoordinatesBig);
            Acquisition::getSettingsPerShot(sourceSettings, sourceSettingsBig, cutCoordinates);
        } else {
            sourceSettings = sources.getSourceSettings(); 
        }
        std::vector<scai::IndexType> uniqueShotNos;
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
        IndexType numshotsIncr = sourceSettingsEncode.size();
        std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncodeTemp; // to control the useSourceEncode
        KITGPI::Acquisition::Receivers<ValueType> receivers;
        
        taperEncode.clear();
        for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
            if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumberEncode) {
                IndexType shotNumber = uniqueShotNos[sourceSettingsEncode[shotInd].row];                 
                std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
                Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
                sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);
                
                if (config.get<IndexType>("useReceiversPerShot") != 0) {
                    receivers.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncodeTemp);
                }
                
                init(dist, ctx, sources, receivers, config, modelCoordinates);  
                taperEncode.push_back(taper);
            }
        }
        end_t_shot = common::Walltime::get();
        HOST_PRINT(commShot, "Shot number " << shotNumberEncode << ": Calculate encoded source and receiver taper in " << end_t_shot - start_t_shot << " sec.\n");
    }
}

template class KITGPI::Preconditioning::SourceReceiverTaper<double>;
template class KITGPI::Preconditioning::SourceReceiverTaper<float>;
