

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
            ZeroLagXcorr3Delastic(){equationType="elastic"; numDimension=3;};

            //! Default destructor
            ~ZeroLagXcorr3Delastic(){};

            explicit ZeroLagXcorr3Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow);

            void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            int getNumDimension() const;
            std::string getEquationType() const;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void write(std::string filename, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

          private:
              
            using ZeroLagXcorr<ValueType>::numDimension;
            using ZeroLagXcorr<ValueType>::equationType; 
            
            /* required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrRho;
            using ZeroLagXcorr<ValueType>::xcorrLambda;
            using ZeroLagXcorr<ValueType>::xcorrMuA;
            using ZeroLagXcorr<ValueType>::xcorrMuB;
            using ZeroLagXcorr<ValueType>::xcorrMuC;

            std::string type;
        };
    }
}
