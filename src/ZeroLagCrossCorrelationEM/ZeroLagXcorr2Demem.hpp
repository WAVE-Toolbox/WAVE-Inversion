

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

        /*! \brief Class to hold and calculate the zero lag cross correlated wavefields for 2Dememgradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Demem : public ZeroLagXcorrEM<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Demem(){equationType="emem"; numDimension=2;};

            //! Default destructor
            ~ZeroLagXcorr2Demem(){};

            explicit ZeroLagXcorr2Demem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config);

            void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void gatherWavefields(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep) override;
            void sumWavefields(KITGPI::Workflow::Workflow<ValueType> const &workflow, ValueType DT) override;
            void applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
            int getNumDimension() const;
            std::string getEquationType() const;
            
            /* Getter routines for non-required wavefields: Will throw an error */    
            scai::lama::DenseVector<ValueType> const &getXcorrREpsilonSigmaEM() const override;   

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config) override;

            void write(std::string filename, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
          private:
              
            using ZeroLagXcorr<ValueType>::numDimension;
            using ZeroLagXcorr<ValueType>::equationType;
            using ZeroLagXcorr<ValueType>::gradientKernel;
            using ZeroLagXcorr<ValueType>::gradientDomain;
            using ZeroLagXcorr<ValueType>::decomposition;
            using ZeroLagXcorr<ValueType>::dtinversion;
            using ZeroLagXcorr<ValueType>::dhinversion;
             
            
            /* required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrSigmaEM;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEM;
            
            /* non required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrREpsilonSigmaEM;
            
            std::string type;
        };
    }
}
