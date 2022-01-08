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
            
            IndexType NR = receiverIndices.size();
            lama::DenseVector<ValueType> offsetSign(NR, 0.0);
            if(recY[NR-1] != recY[0] && abs(recX[NR-1] - recX[0]) / abs(recY[NR-1] - recY[0]) < 1){ // borehole
                offsetSign = recY - sourceCoords.y;
            } else { // surface
                offsetSign = recX - sourceCoords.x;
            }
            offsetSign.unaryOp(offsetSign, common::UnaryOp::SIGN);
            offsetSign += 0.1; // except for the same depth source and receiver in borehole case.
            offsetSign.unaryOp(offsetSign, common::UnaryOp::SIGN);
            result *= offsetSign;
        }  
        
        /*! \brief calculate Gaussian window
        \param vector2D gradient vector
        \param GaussianKernel Gaussian kernel
        \param PX size of kernel
        \param modelCoordinates coordinate class object of the model
        \param velocityMean mean velocity of the model
        \param FCmax max frequency
        */
        template <typename ValueType>
        void calcGaussianKernelFor2DVector(scai::lama::DenseVector<ValueType> const vector2D, scai::lama::CSRSparseMatrix<ValueType> &GaussianKernel, IndexType &PX, IndexType &PY, IndexType NX, IndexType NY, ValueType DH, ValueType velocityMean, ValueType FCmax, IndexType smoothGradient)
        {
            scai::hmemo::ContextPtr ctx = vector2D.getContextPtr();
			
			/* define filter size as fraction of reference velocity wavelength */
            ValueType wavelengthMin = velocityMean / FCmax;
            ValueType wavelengthFraction = smoothGradient % 10;
            smoothGradient = (smoothGradient - wavelengthFraction) / 10;
            ValueType wavelengthFractionX = 1.0 / wavelengthFraction;
            ValueType wavelengthFractionY = 0;
            if (smoothGradient == 2) {
                wavelengthFractionY = 1.0 / wavelengthFraction;
            }
            PX = round(wavelengthFractionX * wavelengthMin / DH);
            PY = round(wavelengthFractionY * wavelengthMin / DH);
			
            if (!(PX % 2)) {
                PX += 1;
            }
            if (!(PY % 2)) {
                PY += 1;
            }
    
            IndexType PXhalf = PX / 2;
            IndexType PYhalf = PY / 2;
                    
			ValueType sigmaX = PXhalf / 2.0;
			ValueType sigmaX2 = 2.0 * sigmaX * sigmaX;
			ValueType sigmaY = PYhalf / 2.0;
			ValueType sigmaY2 = 2.0 * sigmaY * sigmaY;
            if (sigmaX2 != 0)
                sigmaX2 = 1.0 / sigmaX2;
            if (sigmaY2 != 0)
                sigmaY2 = 1.0 / sigmaY2;
            
			/* create filter kernel */
            scai::lama::DenseVector<ValueType> kernel(PX*PY, 0, ctx);
			for (IndexType ix=-PXhalf; ix<=PXhalf; ix++){
			    for (IndexType iy=-PYhalf; iy<=PYhalf; iy++){						      
			        kernel[(iy+PYhalf)*PX+ix+PXhalf] = exp(-((ix*ix)*sigmaX2)-((iy*iy)*sigmaY2));      
			    }
			}
            kernel *= 1 / kernel.sum();
            
			/* create filter matrix */
            dmemo::DistributionPtr dist = vector2D.getDistributionPtr();
            dmemo::DistributionPtr distPadded(new scai::dmemo::NoDistribution((NX+PX-1)*(NY+PY-1)));
            GaussianKernel.allocate(dist, distPadded);
            GaussianKernel.setContextPtr(ctx);
            hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
            IndexType colIndex;
            lama::MatrixAssembly<ValueType> assembly;
            dist->getOwnedIndexes(ownedIndexes);
            IndexType X;
            IndexType Y;
            for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
                Y = ownedIndex / NX;
                X = ownedIndex - NX * Y;
                for (IndexType ix=-PXhalf; ix<=PXhalf; ix++){
                    for (IndexType iy=-PYhalf; iy<=PYhalf; iy++){					 
                        colIndex=(Y+PYhalf+iy)*(NX+PX-1)+X+PXhalf+ix;
                        assembly.push(ownedIndex, colIndex, kernel[(iy+PYhalf)*PX+ix+PXhalf]);
                    }
                }
            }
            GaussianKernel.fillFromAssembly(assembly);
        }
        
        /*! \brief pad a 2D vector
        \param vector2D gradient vector
        \param PX size of kernel
        \param modelCoordinates coordinate class object of the model
        */
        template <typename ValueType>
        void pad2DVector(scai::lama::DenseVector<ValueType> const vector2D, scai::lama::DenseVector<ValueType> &vector2Dpadded, IndexType NX, IndexType NY, IndexType PX, IndexType PY)
        {
            scai::hmemo::ContextPtr ctx = vector2D.getContextPtr();
			
            if (!(PX % 2)) {
                PX += 1;
            }
            if (!(PY % 2)) {
                PY += 1;
            }
            
            IndexType PXhalf = PX / 2;
            IndexType PYhalf = PY / 2;
            scai::lama::DenseVector<ValueType> vector2DTemp((NX+PX-1)*(NY+PY-1), 0, ctx);
            vector2Dpadded = vector2DTemp;
                        
			/* center */
            for (IndexType iy = 0; iy < NY; iy++) {
                for (IndexType ix = 0; ix < NX; ix++) {
                    vector2Dpadded[(PYhalf+iy)*(NX+PX-1)+PXhalf+ix]=vector2D[iy*NX+ix];
                }
            }
			/* left and right */
            for (IndexType iy = 0; iy < NY; iy++) {
                for (IndexType ix = 0; ix < PXhalf; ix++) {
                    vector2Dpadded[(PYhalf+iy)*(NX+PX-1)+ix]=vector2D[iy*NX];
                    vector2Dpadded[(PYhalf+iy)*(NX+PX-1)+NX+PXhalf+ix]=vector2D[iy*NX+NX-1];
                }
            }
			/* top and bottom */
            for (IndexType iy = 0; iy < PYhalf; iy++) {
                for (IndexType ix = 0; ix < NX+PX-1; ix++) {
                    vector2Dpadded[iy*(NX+PX-1)+ix]=vector2Dpadded[PYhalf*(NX+PX-1)+ix];
                    vector2Dpadded[(NY+PYhalf+iy)*(NX+PX-1)+ix]=vector2Dpadded[(NY+PYhalf-1)*(NX+PX-1)+ix];
                }
            }
        }
    }    
}
