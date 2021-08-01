

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

        /*! \brief Class to hold and calculate the zero lag cross correlated wavefields for 2Dviscotmemgradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Dviscotmem : public ZeroLagXcorrEM<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Dviscotmem(){equationType="viscotmem"; numDimension=2;};

            //! Default destructor
            ~ZeroLagXcorr2Dviscotmem(){};

            explicit ZeroLagXcorr2Dviscotmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow);

            void resetXcorr(KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) override;

            void update(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow, scai::IndexType gradientType, scai::IndexType decomposeType) override;
            void updateHessianVectorProduct(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, Wavefields::WavefieldsEM<ValueType> &forwardWavefield2ndOrder, Wavefields::WavefieldsEM<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) override;
            int getNumDimension() const;
            std::string getEquationType() const;
            
            /* Getter routines for non-required wavefields: Will throw an error */     

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) override;

            void write(std::string filename, scai::IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) override;

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
