#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/GridReadAccess.hpp>

#include <Acquisition/AcquisitionGeometry.hpp>
#include <Configuration/Configuration.hpp>

#include "../Gradient/GradientFactory.hpp"

namespace KITGPI
{
    
    //! \brief Gradient preconditioning namespace
    namespace Preconditioning
    {
        
        /*! \brief Class to damp a gradient in the vicinity of the sources and receivers
         * 
         */
        template <typename ValueType>
        class SourceReceiverTaper
        {
        
        public:
            /* Default constructor and destructor */
            SourceReceiverTaper(){};
            ~SourceReceiverTaper(){};

            scai::lama::SparseVector<ValueType> getTaper();
            std::vector<scai::lama::SparseVector<ValueType>> getTaperEncode();
            void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Acquisition::AcquisitionGeometry<ValueType> const &Acquisition, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType radius);            
            void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Acquisition::AcquisitionGeometry<ValueType> const &sources, KITGPI::Acquisition::AcquisitionGeometry<ValueType> const &receivers, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates);      
            void init(scai::dmemo::CommunicatorPtr commShot, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, scai::IndexType shotNumberEncode);
            
            void apply(KITGPI::Gradient::Gradient<ValueType> &gradient);
            void stack(scai::lama::SparseVector<ValueType> taperSingle);
                        
        private:
            scai::lama::SparseVector<ValueType> taper;
            std::vector<scai::lama::SparseVector<ValueType>> taperEncode;
        };
    }
}
