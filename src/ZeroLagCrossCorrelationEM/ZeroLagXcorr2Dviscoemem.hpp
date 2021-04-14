

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "ZeroLagXcorr.hpp"

#include "../WorkflowEM/Workflow.hpp"

namespace KITGPI
{

    namespace ZeroLagXcorr
    {

        /*! \brief Class to hold and calculate the zero lag cross correlated wavefieldsEM for 2Dviscoememgradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Dviscoemem : public ZeroLagXcorrEM<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Dviscoemem(){equationTypeEM="viscoemem"; numDimension=2;};

            //! Default destructor
            ~ZeroLagXcorr2Dviscoemem(){};

            explicit ZeroLagXcorr2Dviscoemem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM);

            void resetXcorr(KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) override;

            void update(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) override;
            void updateHessianVectorProduct(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, Wavefields::WavefieldsEM<ValueType> &forwardWavefield2ndOrder, Wavefields::WavefieldsEM<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) override;
            int getNumDimension() const;
            std::string getEquationType() const;
            
            /* Getter routines for non-required wavefieldsEM: Will throw an error */    

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) override;

            void write(std::string type, scai::IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) override;
            void writeSnapshot(scai::IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM);

          private:
              
            using ZeroLagXcorrEM<ValueType>::numDimension;
            using ZeroLagXcorrEM<ValueType>::equationTypeEM;
            
            /* required wavefieldsEM */
            using ZeroLagXcorrEM<ValueType>::xcorrSigmaEM;
            using ZeroLagXcorrEM<ValueType>::xcorrEpsilonEM;
            using ZeroLagXcorrEM<ValueType>::xcorrRSigmaEM;
            using ZeroLagXcorrEM<ValueType>::xcorrREpsilonEM;
            
            /* non required wavefieldsEM */
            
            std::string type = "ViscoEMEM2D";
        };
    }
}
