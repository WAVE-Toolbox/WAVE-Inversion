

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

        /*! \brief Class to hold and calculate the zero lag cross correlated wavefieldsEM for 2Dememgradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Demem : public ZeroLagXcorrEM<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Demem(){equationTypeEM="emem"; numDimension=2;};

            //! Default destructor
            ~ZeroLagXcorr2Demem(){};

            explicit ZeroLagXcorr2Demem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM);

            void resetXcorr(KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) override;

            void update(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) override;
            void updateHessianVectorProduct(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, Wavefields::WavefieldsEM<ValueType> &forwardWavefield2ndOrder, Wavefields::WavefieldsEM<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) override;
            int getNumDimension() const;
            std::string getEquationType() const;
            
            /* Getter routines for non-required wavefieldsEM: Will throw an error */    
            scai::lama::DenseVector<ValueType> const &getXcorrRSigmaEM() const override;   
            scai::lama::DenseVector<ValueType> const &getXcorrREpsilonEM() const override;  

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
            
            /* non required wavefieldsEM */
            using ZeroLagXcorrEM<ValueType>::xcorrRSigmaEM;
            using ZeroLagXcorrEM<ValueType>::xcorrREpsilonEM;
            
            std::string type = "EMEM2D";
        };
    }
}
