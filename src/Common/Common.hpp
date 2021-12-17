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
        
        /*! \brief calculate Gaussian window
        \param vector2D gradient vector
        \param GaussianKernel Gaussian kernel
        \param ksize size of kernel
        \param modelCoordinates coordinate class object of the model
        \param velocityMean mean velocity of the model
        \param FCmax max frequency
        */
        template <typename ValueType>
        void calcGaussianKernelFor2DVector(scai::lama::DenseVector<ValueType> &vector2D, scai::lama::CSRSparseMatrix<ValueType> &GaussianKernel, IndexType &ksize, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType velocityMean, ValueType FCmax)
        {
            IndexType NX = modelCoordinates.getNX();
            IndexType NY = modelCoordinates.getNY();
            scai::hmemo::ContextPtr ctx = vector2D.getContextPtr();
			
			/* define filter size as fraction of reference velocity wavelength */
            ValueType wavelengthMin = velocityMean / FCmax;
            ValueType WD_DAMP = 0.4;
            ksize = round(WD_DAMP * wavelengthMin / modelCoordinates.getDH());
			
            if (!(ksize % 2)) {
                ksize += 1;
            }
    
            IndexType khalf = ksize / 2;
                    
			ValueType sigma = khalf / 2.0;
			ValueType sigma2 = 2.0 * sigma * sigma;
            
			/* create filter kernel */
            scai::lama::DenseVector<ValueType> kernel(ksize*ksize, 0, ctx);
			for (IndexType ix=-khalf; ix<=khalf; ix++){
			    for (IndexType iy=-khalf; iy<=khalf; iy++){						      
			        kernel[(iy+khalf)*ksize+ix+khalf] = exp(-((ix*ix)/sigma2) - ((iy*iy)/sigma2));      
			    }
			}
            kernel *= 1 / kernel.sum();
            
			/* create filter matrix */
            dmemo::DistributionPtr dist = vector2D.getDistributionPtr();
            dmemo::DistributionPtr distPadded(new scai::dmemo::NoDistribution((NX+ksize-1)*(NY+ksize-1)));
            GaussianKernel.allocate(dist, distPadded);
            GaussianKernel.setContextPtr(ctx);
            hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
            IndexType colIndex;
            lama::MatrixAssembly<ValueType> assembly;
            dist->getOwnedIndexes(ownedIndexes);
            KITGPI::Acquisition::coordinate3D coordinate;
            for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
                coordinate = modelCoordinates.index2coordinate(ownedIndex);
                for (IndexType ix=-khalf; ix<=khalf; ix++){
                    for (IndexType iy=-khalf; iy<=khalf; iy++){					 
                        colIndex=(coordinate.y+khalf+iy)*(NX+ksize-1)+coordinate.x+khalf+ix;
                        assembly.push(ownedIndex, colIndex, kernel[(iy+khalf)*ksize+ix+khalf]);
                    }
                }
            }
            GaussianKernel.fillFromAssembly(assembly);
        }
        
        /*! \brief pad a 2D vector
        \param vector2D gradient vector
        \param ksize size of kernel
        \param modelCoordinates coordinate class object of the model
        */
        template <typename ValueType>
        void pad2DVector(scai::lama::DenseVector<ValueType> &vector2D, scai::lama::DenseVector<ValueType> &vector2Dpadded, IndexType ksize, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates)
        {
            IndexType NX = modelCoordinates.getNX();
            IndexType NY = modelCoordinates.getNY();
            scai::hmemo::ContextPtr ctx = vector2D.getContextPtr();
			
            if (!(ksize % 2)) {
                ksize += 1;
            }
            
            IndexType khalf = ksize / 2;
            scai::lama::DenseVector<ValueType> vector2DTemp((NX+ksize-1)*(NY+ksize-1), 0, ctx);
            vector2Dpadded = vector2DTemp;
                        
			/* center */
            for (IndexType iy = 0; iy < NY; iy++) {
                for (IndexType ix = 0; ix < NX; ix++) {
                    vector2Dpadded[(khalf+iy)*(NX+ksize-1)+ksize+ix]=vector2D[iy*NX+ix];
                }
            }
			/* left and right */
            for (IndexType iy = 0; iy < NY; iy++) {
                for (IndexType ix = 0; ix < khalf; ix++) {
                    vector2Dpadded[(khalf+iy)*(NX+ksize-1)+ix]=vector2D[iy*NX];
                    vector2Dpadded[(khalf+iy)*(NX+ksize-1)+NX+khalf+ix]=vector2D[iy*NX+NX-1];
                }
            }
			/* top and bottom */
            for (IndexType iy = 0; iy < khalf; iy++) {
                for (IndexType ix = 0; ix < NX+ksize-1; ix++) {
                    vector2Dpadded[iy*(NX+ksize-1)+ix]=vector2Dpadded[khalf*(NX+ksize-1)+ix];
                    vector2Dpadded[(NY+khalf+iy)*(NX+ksize-1)+ix]=vector2Dpadded[(NY+khalf-1)*(NX+ksize-1)+ix];
                }
            }
        }
    }    
}
