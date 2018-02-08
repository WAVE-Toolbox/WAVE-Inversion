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

            virtual scai::lama::Scalar calc(KITGPI::Acquisition::Receivers<ValueType> const &receivers1, KITGPI::Acquisition::Receivers<ValueType> const &receivers2) = 0;
            virtual scai::lama::Scalar calc(KITGPI::Acquisition::SeismogramHandler<ValueType> const &seismoHandler1, KITGPI::Acquisition::SeismogramHandler<ValueType> const &seismoHandler2) = 0;
            virtual scai::lama::Scalar calc(KITGPI::Acquisition::Seismogram<ValueType> &seismogram1, KITGPI::Acquisition::Seismogram<ValueType> &seismogram2) = 0;
            scai::lama::Scalar getMisfitSum(int iteration);
            scai::lama::Scalar getMisfitShot(int iteration, int shotNumber);
            void add(scai::lama::DenseVector<ValueType> vector);

            //! \brief Misfit pointer
            typedef std::shared_ptr<Misfit<ValueType>> MisfitPtr;
            
        protected:
            
            /* Default constructor and destructor */
            Misfit(){};
            ~Misfit(){};
            
            std::vector<scai::lama::DenseVector<ValueType>> misfitShot;
            
        };
    }
}
