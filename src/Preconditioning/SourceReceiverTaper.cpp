#include "SourceReceiverTaper.hpp"
using namespace scai;

template <typename ValueType>
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::getTaper()
{
    //return ;
}

template <typename ValueType>
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::apply(KITGPI::Gradient::Gradient<ValueType> &gradient)
{   
    gradient *= taper;
}

template <typename ValueType>
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::apply(KITGPI::Gradient::GradientEM<ValueType> &gradientEM)
{   
    gradientEM *= taper;
}

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

    // distributed acquistion coordinates needed for each processor
    auto acquisition1DCoordinates = Acquisition.get1DCoordinates();
    acquisition1DCoordinates.replicate();
    auto rAcquisition1DCoordinates = scai::hmemo::hostReadAccess(acquisition1DCoordinates.getLocalValues());

    if(radius>0){
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

                KITGPI::Acquisition::coordinate3D centerCoord;

                // get 3D coordinate of source or receiver position
                centerCoord = modelCoordinates.index2coordinate(SourceCoordinate);
                
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
                    }
                    else {
                        value=0;
                    }
                    
                assembly.push(ownedIndex, value);
                }
            }
        }
    }

    taper.setContextPtr(ctx);
    taper.fillFromAssembly(assembly);
}

template <typename ValueType>
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Acquisition::AcquisitionGeometryEM<ValueType> const &Acquisition, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType radius)
{
    /* Get local "global" indices */
    hmemo::HArray<IndexType> ownedIndexes;
    dist->getOwnedIndexes(ownedIndexes); // get global indexes for local part

    //gets coordinate index for first source to be changed
    IndexType SourceCoordinate;
    
    IndexType taperType = config.get<IndexType>("sourceReceiverTaperType");

    taper.setSameValue(dist, 1.0); // sparse vector taper gets zero element 1.0

    lama::VectorAssembly<ValueType> assembly;

    // distributed acquistion coordinates needed for each processor
    auto acquisition1DCoordinates = Acquisition.get1DCoordinates();
    acquisition1DCoordinates.replicate();
    auto rAcquisition1DCoordinates = scai::hmemo::hostReadAccess(acquisition1DCoordinates.getLocalValues());

    if(radius>0){
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

                KITGPI::Acquisition::coordinate3D centerCoord;

                // get 3D coordinate of source or receiver position
                centerCoord = modelCoordinates.index2coordinate(SourceCoordinate);
                
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
                    }
                    else {
                        value=0;
                    }
                    
                assembly.push(ownedIndex, value);
                }
            }
        }
    }

    taper.setContextPtr(ctx);
    taper.fillFromAssembly(assembly);
}

template class KITGPI::Preconditioning::SourceReceiverTaper<double>;
template class KITGPI::Preconditioning::SourceReceiverTaper<float>;
