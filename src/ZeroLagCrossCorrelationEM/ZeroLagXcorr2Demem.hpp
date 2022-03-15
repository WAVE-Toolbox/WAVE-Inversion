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

            explicit ZeroLagXcorr2Demem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot);

            void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void gatherWavefields(Wavefields::Wavefields<ValueType> &wavefields, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint) override;
            void sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sourceFC, ValueType DT, scai::IndexType shotNumber) override;
            void applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
            int getNumDimension() const;
            std::string getEquationType() const;
            
            /* Getter routines for non-required wavefields: Will throw an error */    
            scai::lama::DenseVector<ValueType> const &getXcorrREpsilonSigmaEM() const override;   

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot) override;

            void write(std::string filename, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
          private:

            typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;
    
            using ZeroLagXcorr<ValueType>::numDimension;
            using ZeroLagXcorr<ValueType>::equationType;
            using ZeroLagXcorr<ValueType>::gradientKernel;
            using ZeroLagXcorr<ValueType>::gradientDomain;
            using ZeroLagXcorr<ValueType>::normalizeGradient;
            using ZeroLagXcorr<ValueType>::useSourceEncode;
            using ZeroLagXcorr<ValueType>::decomposition;
            using ZeroLagXcorr<ValueType>::NT;             
            
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
            using ZeroLagXcorr<ValueType>::EXforward;
            using ZeroLagXcorr<ValueType>::EXadjoint;
            using ZeroLagXcorr<ValueType>::EYforward;
            using ZeroLagXcorr<ValueType>::EYadjoint;
            using ZeroLagXcorr<ValueType>::fEXforward;
            using ZeroLagXcorr<ValueType>::fEXadjoint;
            using ZeroLagXcorr<ValueType>::fEYforward;
            using ZeroLagXcorr<ValueType>::fEYadjoint;
            
            /* non required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrREpsilonSigmaEM;
            
            std::string type;
        };
    }
}
