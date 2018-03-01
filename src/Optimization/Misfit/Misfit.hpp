#pragma once

#include <scai/lama.hpp>
#include <vector>
#include <Acquisition/Receivers.hpp>

namespace KITGPI
{
    
    //! \brief Misfit namespace
    namespace Misfit
    {
        //! \brief Abstract class for the misfit
        /*!
         * As this class is an abstract class, all constructors are protected.
         */
        
        template <typename ValueType>
        class Misfit
        {

        public:

            virtual scai::lama::Scalar calc(KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs) = 0;
            virtual scai::lama::Scalar calc(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs) = 0;
            
            virtual void calcAdjointSources(KITGPI::Acquisition::Receivers<ValueType> &adjointSources, KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs) = 0;
            virtual void calcAdjointSeismogram(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs) = 0;
            
            scai::lama::Scalar getMisfitSum(int iteration);
            scai::lama::Scalar getMisfitShot(int iteration, int shotNumber);
            void addToStorage(scai::lama::DenseVector<ValueType> vector);

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
