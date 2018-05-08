

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "ZeroLagXcorr.hpp"

#include "../../Workflow/Workflow.hpp"

namespace KITGPI
{

    namespace ZeroLagXcorr
    {

        /*! \brief Class to hold and caclulate the zero lag cross correlated wavefields for 2D acoustic gradients
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

            explicit ZeroLagXcorr2Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow);

            void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> const &getShearStress() const override;
            scai::lama::DenseVector<ValueType> const &getNormalStressDiff() const override;
            scai::lama::DenseVector<ValueType> const &getNormalStressSum() const override;
            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void write(std::string type, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void writeSnapshot(scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow);

          private:
            /* required wavefields */
            using ZeroLagXcorr<ValueType>::VSum;
            using ZeroLagXcorr<ValueType>::P;
            /* non required wavefields */
            using ZeroLagXcorr<ValueType>::ShearStress;
            using ZeroLagXcorr<ValueType>::NormalStressDiff;
            using ZeroLagXcorr<ValueType>::NormalStressSum;
            std::string type = "Acoustic2D";
        };
    }
}
