#pragma once

#include <scai/lama.hpp>
#include <vector>
#include <Acquisition/Receivers.hpp>
#include <Common/Hilbert.hpp>
#include "../Common/FK.hpp"
#include <scai/lama/fft.hpp>

namespace KITGPI
{
    
    //! \brief Misfit namespace
    namespace Misfit
    {
        /*! \brief Abstract class to store and calculate the misfit
         *
         * As this class is an abstract class, all constructors are protected.
         */        
        template <typename ValueType>
        class Misfit
        {

        public:            
            /* Default constructor and destructor */
            Misfit(){};
            ~Misfit(){};
            //! \brief Misfit pointer
            typedef std::shared_ptr<Misfit<ValueType>> MisfitPtr;

            virtual void init(KITGPI::Configuration::Configuration config, std::vector<scai::IndexType> misfitTypeHistory, scai::IndexType numshots, ValueType fmax, ValueType vmin) = 0;
            virtual ValueType calc(KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs, scai::IndexType shotInd) = 0;
            
            virtual void calcAdjointSources(KITGPI::Acquisition::Receivers<ValueType> &adjointSources, KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs, scai::IndexType shotInd) = 0;            
            
            void calcReflectSources(KITGPI::Acquisition::Receivers<ValueType> &sourcesReflect, scai::lama::DenseVector<ValueType> reflectivity);  
            
            scai::lama::Vector<ValueType> const &getModelDerivativeX();
            scai::lama::Vector<ValueType> const &getModelDerivativeY();
            void setModelDerivativeX(scai::lama::Vector<ValueType> const &setModelDerivativeX);
            void setModelDerivativeY(scai::lama::Vector<ValueType> const &setModelDerivativeY);
            ValueType calcStablizingFunctionalPerModel(scai::lama::DenseVector<ValueType> modelResidualVec, ValueType focusingParameter, int stablizingFunctionalType);
            
            ValueType getMisfitSum(int iteration);
            scai::lama::DenseVector<ValueType> getMisfitIt(int iteration);
            ValueType getMisfitShot(int iteration, int shotInd);
            void addToStorage(scai::lama::DenseVector<ValueType> vector);
            void addToCrossGradientMisfitStorage(ValueType crossGradientMisfit); 
            ValueType getCrossGradientMisfit(int iteration);
            void clearStorage();            
            
            virtual void appendMisfitTypeShotsToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration) = 0;
            virtual void appendMisfitPerShotToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration) = 0;
            virtual void appendMultiMisfitsToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration) = 0;
            virtual void sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot) = 0;
            scai::lama::DenseVector<ValueType> getMisfitTypeShots();
            void setMisfitTypeShots(scai::lama::DenseVector<ValueType> setMisfitTypeShots);
            std::vector<scai::lama::DenseVector<ValueType>> getMisfitSum0Ratio();
            void setMisfitSum0Ratio(std::vector<scai::lama::DenseVector<ValueType>> setMisfitSum0Ratio);
            ValueType getMisfitResidualMax(int iteration1, int iteration2);
            
            std::string misfitType = "l2";
            std::string multiMisfitType = "l2567891";
            bool saveMultiMisfits = false;
            KITGPI::FK<ValueType> fkHandler;
            scai::IndexType nFFT;
            
            KITGPI::Misfit::Misfit<ValueType> &operator=(KITGPI::Misfit::Misfit<ValueType> const &rhs);
                                    
        protected:
            
            std::vector<scai::lama::DenseVector<ValueType>> misfitStorage;            
            std::vector<scai::lama::DenseVector<ValueType>> misfitStorageL2;          
            std::vector<scai::lama::DenseVector<ValueType>> misfitSum0Ratio;
            std::vector<ValueType> crossGradientMisfitStorage; 
            scai::lama::DenseVector<ValueType> misfitTypeShots; 
            scai::IndexType numMisfitTypes = 6;
            scai::IndexType gradientType = 0;
            std::vector<scai::IndexType> uniqueMisfitTypes{2, 5, 6, 7, 8, 9};
            scai::lama::DenseVector<ValueType> modelDerivativeX; //!< Vector storing model derivative in x direction.
            scai::lama::DenseVector<ValueType> modelDerivativeY; //!< Vector storing model derivative in y direction.
            
        };
    }
}
