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
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Acquisition::AcquisitionGeometry<ValueType> const &Acquisition, KITGPI::Configuration::Configuration config,KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates,IndexType radius)
{
    
    lama::DenseVector<ValueType> taperTmp(dist, 0.0);

    /* Get local "global" indices */
    hmemo::HArray<IndexType> globalIndices;
    dist->getOwnedIndexes(globalIndices); // get global indeces for local part

    IndexType numLocalIndices = globalIndices.size(); // Number of local indices

    //temp taper arrays of local model size
    // Harray initialized with 0.0
    hmemo::HArray<ValueType> taperValues(numLocalIndices, 0.0, ctx);
    hmemo::WriteAccess<ValueType> write_taperValues(taperValues);
    ValueType write_taperValues_temp[numLocalIndices];

    //gets coordinate index for first source to be changed
    IndexType SourceCoordinate;
    
    IndexType taperType = config.get<IndexType>("sourceReceiverTaperType");

    
    //loop over all shots or receivers
    for (IndexType srcRecNum = 0; srcRecNum < Acquisition.get1DCoordinates().size(); srcRecNum++) {

	// get 1D coordinate
        SourceCoordinate = Acquisition.get1DCoordinates().getValue(srcRecNum);

        KITGPI::Acquisition::coordinate3D coord;
        KITGPI::Acquisition::coordinate3D centerCoord;

	// get 3D coordinate of source or receiver position
        centerCoord = modelCoordinates.index2coordinate(SourceCoordinate);
        ValueType distance;

        hmemo::ReadAccess<IndexType> read_globalIndices(globalIndices); // Get read access to localy stored global indices

	if(radius>0){
        //loop over all gridpoint -> this is costly maybe there is a better alternative because the size of the taper is known
        for (IndexType i = 0; i < numLocalIndices; i++) {

            //(re)initialize temp with zeros
            write_taperValues_temp[i] = 0;
	    
	    // get 3D coordinate of current local position i
            coord = modelCoordinates.index2coordinate(read_globalIndices[i]);
        
	    // distance between src/rec position and current local gridpoint
            distance = sqrt((centerCoord.x - coord.x) * (centerCoord.x - coord.x) + (centerCoord.y - coord.y) * (centerCoord.y - coord.y) + (centerCoord.z - coord.z) * (centerCoord.z - coord.z));
	    
	    // taper area:
            if (distance <= radius) {
                if (distance >= 1) {
                    //-1 is needed because sparseVector=Densevector sets zero element to 0 for (grad * taper) the zero element should be 1 -> sparsevector+1  
                    switch (taperType) {
                        case 1:
                            //normalized logarithmic function : log(i)/log(max) for log(i)>1 
                            write_taperValues_temp[i] = std::log(distance) / std::log(radius) - 1;
                            break;
                        case 2:
                            //cos^2 function : cos(i*pi/(2*max))^2
                            write_taperValues_temp[i] = common::Math::pow<ValueType>(common::Math::sin<ValueType>(distance * M_PI / (2.0 * radius)), 2.0) - 1;
                    }
                            
                } else {
                    write_taperValues_temp[i] = -1;
                }
                  // update taper values only if there is no smaller tapervalue from an other source or receiver
                if (write_taperValues_temp[i] < write_taperValues[i])
                    write_taperValues[i] = write_taperValues_temp[i];
            }
        }
    }
    }
    write_taperValues.release();
    //taper radius should be input!
    taperTmp.setDenseValues(taperValues);
//     lama::DenseVector<ValueType> taperTmp(dist, std::move(taperValues) ); // problems here during runtime

    // copy assign tmp vector to sparsevector
    taper.setContextPtr(ctx);
    taper = taperTmp;
    // setting zero element to 1 by adding 1 to all values (non zero values were set by value-1)
    taper += 1;
}

template class KITGPI::Preconditioning::SourceReceiverTaper<double>;
template class KITGPI::Preconditioning::SourceReceiverTaper<float>;
