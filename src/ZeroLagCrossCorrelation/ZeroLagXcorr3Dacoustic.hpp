

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

            explicit ZeroLagXcorr3Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow);

            void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> const &getXcorrMuA() const override;
            scai::lama::DenseVector<ValueType> const &getXcorrMuB() const override;
            scai::lama::DenseVector<ValueType> const &getXcorrMuC() const override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void write(std::string type, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void writeSnapshot(scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow);

          private:
            /* required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrRho;
            using ZeroLagXcorr<ValueType>::xcorrLambda;

            /* non-required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrMuA;
            using ZeroLagXcorr<ValueType>::xcorrMuB;
            using ZeroLagXcorr<ValueType>::xcorrMuC;

            std::string type = "Acoustic3D";
        };
    }
}
