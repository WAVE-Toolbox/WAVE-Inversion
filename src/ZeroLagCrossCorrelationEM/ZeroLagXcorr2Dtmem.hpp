

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

        /*! \brief Class to hold and calculate the zero lag cross correlated wavefields for 2Dtmemgradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Dtmem : public ZeroLagXcorrEM<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Dtmem(){equationType="tmem"; numDimension=2;};

            //! Default destructor
            ~ZeroLagXcorr2Dtmem(){};

            explicit ZeroLagXcorr2Dtmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow);

            void resetXcorr(KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) override;
            void updateHessianVectorProduct(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefieldDerivative, Wavefields::Wavefields<ValueType> &adjointWavefield, Wavefields::Wavefields<ValueType> &forwardWavefield2ndOrder, Wavefields::Wavefields<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) override;
            int getNumDimension() const;
            std::string getEquationType() const;
            
            /* Getter routines for non-required wavefields: Will throw an error */ 
            scai::lama::DenseVector<ValueType> const &getXcorrRSigmaEM() const override;   
            scai::lama::DenseVector<ValueType> const &getXcorrREpsilonEM() const override;       

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) override;

            void write(std::string filename, scai::IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) override;
          private:
              
            using ZeroLagXcorrEM<ValueType>::numDimension;
            using ZeroLagXcorrEM<ValueType>::equationType;
            using ZeroLagXcorrEM<ValueType>::gradientType;
            using ZeroLagXcorrEM<ValueType>::decomposeType;
            
            /* required wavefields */
            using ZeroLagXcorrEM<ValueType>::xcorrSigmaEM;
            using ZeroLagXcorrEM<ValueType>::xcorrEpsilonEM;
            using ZeroLagXcorrEM<ValueType>::xcorrSigmaEMSuRu;
            using ZeroLagXcorrEM<ValueType>::xcorrSigmaEMSdRd;
            using ZeroLagXcorrEM<ValueType>::xcorrSigmaEMSuRd;
            using ZeroLagXcorrEM<ValueType>::xcorrSigmaEMSdRu;
            using ZeroLagXcorrEM<ValueType>::xcorrEpsilonEMSuRu;
            using ZeroLagXcorrEM<ValueType>::xcorrEpsilonEMSdRd;
            using ZeroLagXcorrEM<ValueType>::xcorrEpsilonEMSuRd;
            using ZeroLagXcorrEM<ValueType>::xcorrEpsilonEMSdRu;
            using ZeroLagXcorrEM<ValueType>::xcorrSigmaEMstep;
            using ZeroLagXcorrEM<ValueType>::xcorrEpsilonEMstep;
            
            /* non required wavefields */
            using ZeroLagXcorrEM<ValueType>::xcorrRSigmaEM;
            using ZeroLagXcorrEM<ValueType>::xcorrREpsilonEM;
            using ZeroLagXcorrEM<ValueType>::xcorrRSigmaEMstep;
            using ZeroLagXcorrEM<ValueType>::xcorrREpsilonEMstep;
            
            std::string type;
        };
    }
}
