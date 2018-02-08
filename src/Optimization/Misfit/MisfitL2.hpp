#pragma once

#include <scai/lama.hpp>
#include <Acquisition/Acquisition.hpp>
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
            
            void calc();
//             void calc(KITGPI::Acquisition::Seismogram<ValueType> const &seismogram1, KITGPI::Acquisition::Seismogram<ValueType> const &seismogram2);

        private:
            
            using Misfit<ValueType>::misfitShot;
        
        };
        
    }
}
