

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

        /*! \brief The class ZeroLagXcorr2Dacoustic holds and caclulates the zero lag cross correlated wavefields for 2D acoustic gradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Dacoustic : public ZeroLagXcorr<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Dacoustic(){};

            //! Default destructor
            ~ZeroLagXcorr2Dacoustic(){};

            explicit ZeroLagXcorr2Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void reset() override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield) override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void write(std::string type, IndexType t) override;
            void writeSnapshot(IndexType t);

          private:
            using ZeroLagXcorr<ValueType>::invertForVp;
            using ZeroLagXcorr<ValueType>::invertForDensity;
            /* required wavefields */
            using ZeroLagXcorr<ValueType>::VSum;
            using ZeroLagXcorr<ValueType>::P;

            std::string type = "Acoustic2D";
        };
    }
}