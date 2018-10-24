#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <Acquisition/Coordinates.hpp>

#include <cmath>

namespace KITGPI
{
    //! \brief Common namespace
    namespace Common
    {
        
        /*! \brief Calculate the distance from one index coordinate o a set of index coordinates
        \param result Distance vector from each receiver index to the source index
        \param sourceIndex Index coordinate of the source
        \param receiverIndices Index coordinates from the receivers
        \param NX Number of grid points in x-direction
        \param NY Number of grid points in y-direction
        \param NZ Number of grid points in z-direction
        */
        template<typename ValueType>
        void calcOffsets(scai::lama::DenseVector<ValueType> &result, scai::IndexType sourceIndex, scai::lama::DenseVector<scai::IndexType> const &receiverIndices, scai::IndexType NX, scai::IndexType NY, scai::IndexType /*NZ*/) {
            
            Acquisition::Coordinates coords;
            Acquisition::coordinate3D sourceCoords = coords.index2coordinate(sourceIndex, NX, NY, 0);
            
            scai::lama::DenseVector<scai::IndexType> recX;
            scai::lama::DenseVector<scai::IndexType> recY;
            scai::lama::DenseVector<scai::IndexType> recZ;
            scai::lama::DenseVector<scai::IndexType> RecCoordinates;
            
            RecCoordinates = receiverIndices;

            recZ = RecCoordinates / (NX * NY);
            RecCoordinates -= recZ * (NX * NY);

            recY = RecCoordinates / (NX);
            RecCoordinates -= recY * (NX);

            recX = RecCoordinates;
            
            recX -= sourceCoords.x;
            recY -= sourceCoords.y;
            recZ -= sourceCoords.z;
            
            recX.binaryOpScalar(recX, 2, scai::common::BinaryOp::POW, false);
            recY.binaryOpScalar(recY, 2, scai::common::BinaryOp::POW, false);
            recZ.binaryOpScalar(recZ, 2, scai::common::BinaryOp::POW, false);
            
            recX += recY;
            recX += recZ;
            result.unaryOp(scai::lama::eval<scai::lama::DenseVector<ValueType>>(scai::lama::cast<ValueType>(recX)), scai::common::UnaryOp::SQR);
            
        }
    }
    
}