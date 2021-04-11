#pragma once

#include <scai/lama.hpp>
#include <vector>
#include <Acquisition/Receivers.hpp>

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

            virtual ValueType calc(KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs) = 0;
            virtual ValueType calc(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs) = 0;
            
            virtual void calcAdjointSources(KITGPI::Acquisition::Receivers<ValueType> &adjointSources, KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs) = 0;
            virtual void calcAdjointSeismogram(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs) = 0;
            
            ValueType calcStablizingFunctionalPerModel(scai::lama::DenseVector<ValueType> modelResidualVec, ValueType focusingParameter, int stablizingFunctionalType);
            
            ValueType getMisfitSum(int iteration);
            scai::lama::DenseVector<ValueType> getMisfitIt(int iteration);
            ValueType getMisfitShot(int iteration, int shotNumber);
            void addToStorage(scai::lama::DenseVector<ValueType> vector);
            void clearStorage();
            
            //! \brief Misfit pointer
            typedef std::shared_ptr<Misfit<ValueType>> MisfitPtr;
            
        protected:
            
            /* Default constructor and destructor */
            Misfit(){};
            ~Misfit(){};
            
            std::vector<scai::lama::DenseVector<ValueType>> misfitStorage;
            
        };
    }
}
