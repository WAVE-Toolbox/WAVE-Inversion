#pragma once

#include <scai/lama.hpp>
#include <Acquisition/Acquisition.hpp>
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Seismogram.hpp>
#include <AcquisitionEM/Acquisition.hpp>
#include <AcquisitionEM/Receivers.hpp>
#include <AcquisitionEM/Seismogram.hpp>
#include "Misfit.hpp"

namespace KITGPI
{
    
    //! \brief Misfit namespace
    namespace Misfit
    {
        //! \brief Class for least-squares misfit

        template <typename ValueType>
        class MisfitL2 : public Misfit<ValueType>
        {

        public:
            /* Default constructor and destructor */
            MisfitL2(){}; // set default misfitType for step length search.
            ~MisfitL2(){};
            
            void init(KITGPI::Configuration::Configuration config, std::vector<scai::IndexType> misfitTypeHistory, scai::IndexType numshots);
            void appendMisfitTypeShotsToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration);
            void appendMisfitsToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration);
            void sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot);
                   
            ValueType calc(KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs, scai::IndexType shotInd);
            ValueType calc(KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversSyn, KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversObs, scai::IndexType shotInd);            
            void calcAdjointSources(KITGPI::Acquisition::Receivers<ValueType> &adjointSources, KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs, scai::IndexType shotInd);
            void calcAdjointSources(KITGPI::Acquisition::ReceiversEM<ValueType> &adjointSources, KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversSyn, KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversObs, scai::IndexType shotInd);
            
            ValueType calcL2(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);
            ValueType calcL2(KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs);            
            void calcAdjointSeismogramL2(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);            
            void calcAdjointSeismogramL2(KITGPI::Acquisition::SeismogramEM<ValueType> &seismogramAdj, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs);            
            
            ValueType calcL2Normalized(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);
            ValueType calcL2Normalized(KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs);            
            void calcAdjointSeismogramL2Normalized(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);            
            void calcAdjointSeismogramL2Normalized(KITGPI::Acquisition::SeismogramEM<ValueType> &seismogramAdj, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs);
                        
            ValueType calcL2Envelope(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);
            ValueType calcL2Envelope(KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs);            
            void calcAdjointSeismogramL2Envelope(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);            
            void calcAdjointSeismogramL2Envelope(KITGPI::Acquisition::SeismogramEM<ValueType> &seismogramAdj, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs);
                              
        private:            
            using Misfit<ValueType>::misfitStorage;  
            using Misfit<ValueType>::misfitStorageL2;  
            using Misfit<ValueType>::misfitSum0Ratio;   
            using Misfit<ValueType>::misfitType;     
            using Misfit<ValueType>::saveMultiMisfits;     
            using Misfit<ValueType>::misfitTypeShots; 
            using Misfit<ValueType>::uniqueMisfitTypes;
            scai::IndexType numMisfitTypes = 3;
        };        
    }
}
