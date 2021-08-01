#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <Acquisition/Coordinates.hpp>
#include <Common/Common.hpp>

#include <cmath>

using namespace scai;
namespace KITGPI
{
    //! \brief Common namespace
    namespace Common
    {    
        /*! \brief Calculate the distance from one index coordinate o a set of index coordinates
        \param result Distance vector from each receiver index to the source index
        \param sourceIndex Index coordinate of the source
        \param receiverIndices Index coordinates from the receivers
        \param modelCoordinates coordinate class object of the model
        */
        template<typename ValueType>
        void calcOffsets(lama::DenseVector<ValueType> &result, IndexType sourceIndex, lama::DenseVector<IndexType> const &receiverIndices, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates) 
        {
            Acquisition::coordinate3D sourceCoords = modelCoordinates.index2coordinate(sourceIndex);
            
            lama::DenseVector<ValueType> recX;
            lama::DenseVector<ValueType> recY;
            lama::DenseVector<ValueType> recZ;
            lama::DenseVector<ValueType> RecCoordinates;
            IndexType NX = modelCoordinates.getNX();
            IndexType NZ = modelCoordinates.getNZ();
            
            RecCoordinates = lama::cast<ValueType>(receiverIndices);

            recY = lama::floor<ValueType>(lama::eval<lama::DenseVector<ValueType>>(RecCoordinates / (NX * NZ)));
            RecCoordinates -= recY * (NX * NZ);

            recZ = lama::floor<ValueType>(lama::eval<lama::DenseVector<ValueType>>(RecCoordinates / (NX)));
            RecCoordinates -= recZ * (NX);

            recX = RecCoordinates;
            
            recX -= ValueType(sourceCoords.x);
            recY -= ValueType(sourceCoords.y);
            recZ -= ValueType(sourceCoords.z);
            
            recX.binaryOpScalar(recX, 2.0, common::BinaryOp::POW, false);
            recY.binaryOpScalar(recY, 2.0, common::BinaryOp::POW, false);
            recZ.binaryOpScalar(recZ, 2.0, common::BinaryOp::POW, false);
            
            recX += recY;
            recX += recZ;
            result.unaryOp(recX, common::UnaryOp::SQRT);   
            result *= modelCoordinates.getDH();
        }               
    }    
}
