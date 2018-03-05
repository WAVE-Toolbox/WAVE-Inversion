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
            
            scai::lama::Scalar calc(KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs);
            scai::lama::Scalar calc(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);
            
            void calcAdjointSources(KITGPI::Acquisition::Receivers<ValueType> &adjointSources, KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs);
            void calcAdjointSeismogram(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);

        private:
            
            using Misfit<ValueType>::misfitStorage;
        
        };
        
    }
}
