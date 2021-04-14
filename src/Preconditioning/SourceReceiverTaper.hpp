#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/GridReadAccess.hpp>

#include <Acquisition/AcquisitionGeometry.hpp>
#include <AcquisitionEM/AcquisitionGeometry.hpp>
#include <Configuration/Configuration.hpp>

#include "../Gradient/GradientFactory.hpp"
#include "../GradientEM/GradientFactory.hpp"

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

            void getTaper();
            void apply(KITGPI::Gradient::Gradient<ValueType> &gradient);
            void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Acquisition::AcquisitionGeometry<ValueType> const &Acquisition, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType radius);
            
            void apply(KITGPI::Gradient::GradientEM<ValueType> &gradientEM);
            void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Acquisition::AcquisitionGeometryEM<ValueType> const &Acquisition, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType radius);

        private:
            scai::lama::SparseVector<ValueType> taper;
        };
    }
}
