

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "ZeroLagXcorr.hpp"

#include "../Workflow/Workflow.hpp"

namespace KITGPI
{

    namespace ZeroLagXcorr
    {

        /*! \brief Class to hold and caclulate the zero lag cross correlated wavefields for 3D elastic gradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr3Delastic : public ZeroLagXcorr<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr3Delastic(){};

            //! Default destructor
            ~ZeroLagXcorr3Delastic(){};

            explicit ZeroLagXcorr3Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow);

            void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &/*forwardWavefield*/, Wavefields::Wavefields<ValueType> &/*adjointWavefield*/, KITGPI::Workflow::Workflow<ValueType> const &/*workflow*/) override{COMMON_THROWEXCEPTION("elastic convolution is not implemented yet.")};

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> const &getP() const override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void write(std::string type, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void writeSnapshot(scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow);
	    
          private:
            /* required wavefields */
            using ZeroLagXcorr<ValueType>::VSum;

            /* non-required wavefields */
            using ZeroLagXcorr<ValueType>::P;

            std::string type = "Elastic3D";
        };
    }
}
