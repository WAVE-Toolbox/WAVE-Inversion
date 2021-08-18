

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

        /*! \brief Class to hold and calculate the zero lag cross correlated wavefields for 2Dviscoememgradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Dviscoemem : public ZeroLagXcorrEM<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Dviscoemem(){equationType="viscoemem"; numDimension=2;};

            //! Default destructor
            ~ZeroLagXcorr2Dviscoemem(){};

            explicit ZeroLagXcorr2Dviscoemem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow);

            void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void updateHessianVectorProduct(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefieldDerivative, Wavefields::Wavefields<ValueType> &adjointWavefield, Wavefields::Wavefields<ValueType> &forwardWavefield2ndOrder, Wavefields::Wavefields<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            int getNumDimension() const;
            std::string getEquationType() const;
            
            /* Getter routines for non-required wavefields: Will throw an error */    

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void write(std::string filename, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
          private:
              
            using ZeroLagXcorrEM<ValueType>::numDimension;
            using ZeroLagXcorrEM<ValueType>::equationType;
            
            /* required wavefields */
            using ZeroLagXcorrEM<ValueType>::xcorrSigmaEM;
            using ZeroLagXcorrEM<ValueType>::xcorrEpsilonEM;
            using ZeroLagXcorrEM<ValueType>::xcorrRSigmaEM;
            using ZeroLagXcorrEM<ValueType>::xcorrREpsilonEM;
            
            /* non required wavefields */
            
            std::string type;
        };
    }
}
