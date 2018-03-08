

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "ZeroLagXcorr.hpp"

namespace KITGPI
{

    namespace ZeroLagXcorr
    {

        /*! \brief Class to hold and caclulate the zero lag cross correlated wavefields for 3D acoustic gradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr3Dacoustic : public ZeroLagXcorr<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr3Dacoustic(){};

            //! Default destructor
            ~ZeroLagXcorr3Dacoustic(){};

            explicit ZeroLagXcorr3Dacoustic(Configuration::Configuration const &config,scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void resetXcorr() override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield) override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(Configuration::Configuration const &config,scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void write(std::string type, IndexType t) override;
            void writeSnapshot(IndexType t);
	    
	     using ZeroLagXcorr<ValueType>::invertForVp;
            using ZeroLagXcorr<ValueType>::invertForDensity;

          private:
	    
            /* required wavefields */
            using ZeroLagXcorr<ValueType>::VSum;
            using ZeroLagXcorr<ValueType>::P;

            /* non-required wavefields */

            std::string type = "Acoustic3D";
        };
    }
}
