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
    
    //! \brief Preconditioning namespace
    namespace Preconditioning
    {
        template <typename ValueType>
        class SourceReceiverTaper
        {
        
        public:
            /* Default constructor and destructor */
            SourceReceiverTaper(){};
            ~SourceReceiverTaper(){};

            void getTaper();
            void apply(KITGPI::Gradient::Gradient<ValueType> &gradient);
            void init(scai::dmemo::DistributionPtr dist,scai::hmemo::ContextPtr ctx,KITGPI::Acquisition::AcquisitionGeometry<ValueType> const &Acquisition,KITGPI::Configuration::Configuration config,IndexType radius);


        private:
            scai::lama::SparseVector<ValueType> taper;
        };
    }
}
