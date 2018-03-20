

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

        /*! \brief Class to hold and caclulate the zero lag cross correlated wavefields for 2D elastic gradients 
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Delastic : public ZeroLagXcorr<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Delastic(){};

            //! Default destructor
            ~ZeroLagXcorr2Delastic(){};

            explicit ZeroLagXcorr2Delastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void resetXcorr() override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield) override;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> const &getP() const override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void write(std::string type, IndexType t) override;
            void writeSnapshot(IndexType t);

            using ZeroLagXcorr<ValueType>::invertForVp;
            using ZeroLagXcorr<ValueType>::invertForVs;
            using ZeroLagXcorr<ValueType>::invertForDensity;

          private:
            /* required wavefields */
            using ZeroLagXcorr<ValueType>::VSum;
            using ZeroLagXcorr<ValueType>::ShearStress;
            using ZeroLagXcorr<ValueType>::NormalStressDiff;
            using ZeroLagXcorr<ValueType>::NormalStressSum;
            /* non-required wavefields */
            using ZeroLagXcorr<ValueType>::P;

            std::string type = "Elastic2D";
        };
    }
}
