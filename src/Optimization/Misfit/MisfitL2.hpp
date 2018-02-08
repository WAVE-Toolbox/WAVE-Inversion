#pragma once

#include <scai/lama.hpp>
#include <Acquisition/Acquisition.hpp>
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Seismogram.hpp>
#include "Misfit.hpp"


namespace KITGPI
{
    
    //! \brief Misfit namespace
    namespace Misfit
    {
        //! \brief Least-squares misfit class

        template <typename ValueType>
        class MisfitL2 : public Misfit<ValueType>
        {

        public:
            /* Default constructor and destructor */
            MisfitL2(){};
            ~MisfitL2(){};
            
            scai::lama::Scalar calc(KITGPI::Acquisition::Receivers<ValueType> const &receivers1, KITGPI::Acquisition::Receivers<ValueType> const &receivers2);
            scai::lama::Scalar calc(KITGPI::Acquisition::SeismogramHandler<ValueType> const &seismoHandler1, KITGPI::Acquisition::SeismogramHandler<ValueType> const &seismoHandler2);
            scai::lama::Scalar calc(KITGPI::Acquisition::Seismogram<ValueType> &seismogram1, KITGPI::Acquisition::Seismogram<ValueType> &seismogram2);

        private:
            
            using Misfit<ValueType>::misfitShot;
        
        };
        
    }
}
