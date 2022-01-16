

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "ZeroLagXcorrEM.hpp"

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

            explicit ZeroLagXcorr2Dtmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow);

            void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
            int getNumDimension() const;
            std::string getEquationType() const;
            
            /* Getter routines for non-required wavefields: Will throw an error */ 
            scai::lama::DenseVector<ValueType> const &getXcorrREpsilonSigmaEM() const override;     

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void write(std::string filename, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
          private:
              
            using ZeroLagXcorr<ValueType>::numDimension;
            using ZeroLagXcorr<ValueType>::equationType;
            using ZeroLagXcorr<ValueType>::gradientType;
            using ZeroLagXcorr<ValueType>::decomposition;
            
            /* required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrSigmaEM;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEM;
            using ZeroLagXcorr<ValueType>::xcorrSigmaEMSuRu;
            using ZeroLagXcorr<ValueType>::xcorrSigmaEMSdRd;
            using ZeroLagXcorr<ValueType>::xcorrSigmaEMSuRd;
            using ZeroLagXcorr<ValueType>::xcorrSigmaEMSdRu;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEMSuRu;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEMSdRd;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEMSuRd;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEMSdRu;
            using ZeroLagXcorr<ValueType>::xcorrSigmaEMstep;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEMstep;
            
            /* non required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrREpsilonSigmaEM;
            using ZeroLagXcorr<ValueType>::xcorrREpsilonSigmaEMstep;
            
            std::string type;
        };
    }
}
