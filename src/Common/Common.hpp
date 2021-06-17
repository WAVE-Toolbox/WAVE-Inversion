#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <Acquisition/Coordinates.hpp>

#include <cmath>

//   medianfilter.hpp - declarations for 
//   1D and 2D median filter routines
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

#ifndef _MEDIANFILTER_H_
#define _MEDIANFILTER_H_

//   Signal/image element type
typedef float element;

//   1D MEDIAN FILTER, window size 5
//     signal - input signal
//     result - output signal, NULL for inplace processing
//     NX      - length of the signal
void medianfilter(element* signal, element* result, int NX);

//   2D MEDIAN FILTER, window size 3x3
//     image  - input image
//     result - output image, NULL for inplace processing
//     NX      - width of the image
//     NY      - height of the image
void medianfilter(element* image, element* result, int NX, int NY, int ksize, int fillType);

#endif

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
        void calcOffsets(scai::lama::DenseVector<ValueType> &result, scai::IndexType sourceIndex, scai::lama::DenseVector<scai::IndexType> const &receiverIndices, scai::IndexType NX, scai::IndexType NY, scai::IndexType /*NZ*/, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates) 
        {
            Acquisition::coordinate3D sourceCoords = modelCoordinates.index2coordinate(sourceIndex);
            
            scai::lama::DenseVector<ValueType> recX;
            scai::lama::DenseVector<ValueType> recY;
            scai::lama::DenseVector<ValueType> recZ;
            scai::lama::DenseVector<ValueType> RecCoordinates;
            
            RecCoordinates = scai::lama::cast<ValueType>(receiverIndices);

            recZ = scai::lama::floor<ValueType>(scai::lama::eval<scai::lama::DenseVector<ValueType>>(RecCoordinates / (NX * NY)));
            RecCoordinates -= recZ * (NX * NY);

            recY = scai::lama::floor<ValueType>(scai::lama::eval<scai::lama::DenseVector<ValueType>>(RecCoordinates / (NX)));
            RecCoordinates -= recY * (NX);

            recX = RecCoordinates;
            
            recX -= ValueType(sourceCoords.x);
            recY -= ValueType(sourceCoords.y);
            recZ -= ValueType(sourceCoords.z);
            
            recX.binaryOpScalar(recX, 2.0, scai::common::BinaryOp::POW, false);
            recY.binaryOpScalar(recY, 2.0, scai::common::BinaryOp::POW, false);
            recZ.binaryOpScalar(recZ, 2.0, scai::common::BinaryOp::POW, false);
            
            recX += recY;
            recX += recZ;
            result.unaryOp(recX, scai::common::UnaryOp::SQRT);
            
        }
        
        /*! \brief Apply a 2D median filter to model/gradient parameter to filter out the extreme values
        \param vecter2D model/gradient parameter vector
        \param NX Number of grid points in x-direction
        \param NY Number of grid points in y-direction
        \param spatialFDorder spatial FD order used in Derivatives calculation
        */     
        template<typename ValueType>
        void applyMedianFilterTo2DVector(scai::lama::DenseVector<ValueType> &vecter2D, scai::IndexType NX, scai::IndexType NY, scai::IndexType spatialFDorder)
        {    
            element signal_2D[ NY * NX ] = {0};
            for (int i = 0; i < NY * NX; i++) {
                signal_2D[i] = vecter2D.getValue(i);
            }
            element *result;
            result = (element *)malloc(NY * NX * sizeof(element));

            medianfilter(signal_2D, result, NX, NY, spatialFDorder*2+1, 0);
            
            for (int i = 0; i < NY * NX; i++) {
                vecter2D.setValue(i, *(result + i));
            }
        };
    }    
}
