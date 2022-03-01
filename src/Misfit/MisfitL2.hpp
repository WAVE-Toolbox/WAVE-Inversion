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
        //! \brief Class for least-squares misfit

        template <typename ValueType>
        class MisfitL2 : public Misfit<ValueType>
        {

        public:
            /* Default constructor and destructor */
            MisfitL2(){}; // set default misfitType for step length search.
            ~MisfitL2(){};
            
            void init(KITGPI::Configuration::Configuration config, std::vector<scai::IndexType> misfitTypeHistory, scai::IndexType numshots, scai::IndexType useRTM_in, ValueType vmin_in, scai::IndexType &seedtime);
            void appendMisfitTypeShotsToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration);
            void appendMisfitPerShotToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration);
            void appendMultiMisfitsToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration);
            void sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot);
                   
            ValueType calc(KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs, scai::IndexType shotInd) override;          
            void calcAdjointSources(KITGPI::Acquisition::Receivers<ValueType> &adjointSources, KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs, scai::IndexType shotInd) override;
            
            void calcMisfitAndAdjointSources(scai::dmemo::CommunicatorPtr commShot, scai::lama::DenseVector<ValueType> &misfitPerIt, KITGPI::Acquisition::Receivers<ValueType> &adjointSourcesEncode, KITGPI::Acquisition::Receivers<ValueType> const &receivers, KITGPI::Acquisition::Receivers<ValueType> const &receiversTrue, scai::IndexType shotIndTrue, KITGPI::Configuration::Configuration const &config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, ValueType vmin, scai::IndexType &seedtime) override;
            
            ValueType calcL2(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);        
            void calcAdjointSeismogramL2(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);            
            
            ValueType calcL2Convolved(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);        
            void calcAdjointSeismogramL2Convolved(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);    
            
            ValueType calcL2FK(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);        
            void calcAdjointSeismogramL2FK(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);   
            
            ValueType calcL2EnvelopeWeighted(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);  
            void calcAdjointSeismogramL2EnvelopeWeighted(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);            
            
            ValueType calcL2AGC(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);          
            void calcAdjointSeismogramL2AGC(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);            
            
            ValueType calcL2Normalized(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);        
            void calcAdjointSeismogramL2Normalized(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);            
            
            ValueType calcL2Envelope(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);        
            void calcAdjointSeismogramL2Envelope(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);            
                        
            ValueType calcL2InstantaneousPhase(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);        
            void calcAdjointSeismogramL2InstantaneousPhase(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs);            
                              
            using Misfit<ValueType>::misfitType;   
            using Misfit<ValueType>::multiMisfitType;     
            using Misfit<ValueType>::saveMultiMisfits;   
            using Misfit<ValueType>::fkHandler; 
            using Misfit<ValueType>::nFFT; 
            typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;
            ValueType waterLevel = 1e-3;
            
        private:            
            using Misfit<ValueType>::misfitStorage;  
            using Misfit<ValueType>::misfitStorageL2;  
            using Misfit<ValueType>::misfitSum0Ratio;     
            using Misfit<ValueType>::misfitTypeShots; 
            using Misfit<ValueType>::uniqueMisfitTypes;
            using Misfit<ValueType>::numMisfitTypes;
            using Misfit<ValueType>::useRTM;
            using Misfit<ValueType>::writeAdjointSource;
        };        
    }
}
