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
void KITGPI::Preconditioning::SourceReceiverTaper<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Acquisition::AcquisitionGeometry<ValueType> const &Acquisition, KITGPI::Configuration::Configuration config,IndexType radius)
{
    //taper radius should be input!
    lama::DenseVector<ValueType> taperTmp(dist);
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
    lama::Scalar temp;
    IndexType SourceCoordinate;

    
    //loop over all shots or receivers
    for (IndexType srcRecNum = 0; srcRecNum < Acquisition.getCoordinates().size(); srcRecNum++) {

	// get 1D coordinate
        temp = Acquisition.getCoordinates().getValue(srcRecNum);
        SourceCoordinate = temp.getValue<IndexType>();

        KITGPI::Acquisition::Coordinates<ValueType> coordTransform;
        KITGPI::Acquisition::coordinate3D coord;
        KITGPI::Acquisition::coordinate3D centerCoord;

	// get 3D coordinate of source or receiver position
        centerCoord = coordTransform.index2coordinate(SourceCoordinate, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"));
        ValueType distance;

        hmemo::ReadAccess<IndexType> read_globalIndices(globalIndices); // Get read access to localy stored global indices

	if(radius>0){
        //loop over all gridpoint -> this is costly maybe there is a better alternative because the size of the taper is known
        for (IndexType i = 0; i < numLocalIndices; i++) {

            //(re)initialize temp with zeros
            write_taperValues_temp[i] = 0;
	    
	    // get 3D coordinate of current local position i
            coord = coordTransform.index2coordinate(read_globalIndices[i], config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"));
        
	    // distance between src/rec position and current local gridpoint
            distance = sqrt((centerCoord.x - coord.x) * (centerCoord.x - coord.x) + (centerCoord.y - coord.y) * (centerCoord.y - coord.y) + (centerCoord.z - coord.z) * (centerCoord.z - coord.z));
	    
	    // taper area:
            if (distance <= radius) {
                if (distance >= 1) {
                    /*normalized logarithmic function : log(i)/log(max) for log(i)>1 
			-1 is needed because sparseVector=Densevector sets zero element to 0 for (grad * taper) the zero element should be 1 -> sparsevector+1  */
                    write_taperValues_temp[i] = std::log(distance) / std::log(radius) - 1;
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
    taperTmp.setDenseValues(taperValues);

    // copy assign tmp vector to sparsevector
    taper.setContextPtr(ctx);
    taper.allocate(dist);
    taper = taperTmp;
    // setting zero element to 1 by adding 1 to all values (non zero values were set by value-1)
    taper += 1;
}

template class KITGPI::Preconditioning::SourceReceiverTaper<double>;
template class KITGPI::Preconditioning::SourceReceiverTaper<float>;
