

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "ZeroLagXcorrSeismic.hpp"

namespace KITGPI
{

    namespace ZeroLagXcorr
    {

        /*! \brief Class to hold and caclulate the zero lag cross correlated wavefields for 2D acoustic gradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Dacoustic : public ZeroLagXcorrSeismic<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Dacoustic(){equationType="acoustic"; numDimension=2;};

            //! Default destructor
            ~ZeroLagXcorr2Dacoustic(){};

            explicit ZeroLagXcorr2Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot);

            void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void gatherWavefields(Wavefields::Wavefields<ValueType> &wavefields, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint) override;
            void sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sourceFC, ValueType DT, scai::IndexType shotNumber, std::vector<scai::lama::SparseVector<ValueType>> taperEncode) override;
            void applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            int getNumDimension() const;
            std::string getEquationType() const;
            
            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> const &getXcorrMuA() const override;
            scai::lama::DenseVector<ValueType> const &getXcorrMuB() const override;
            scai::lama::DenseVector<ValueType> const &getXcorrMuC() const override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot) override;

            void write(std::string filename, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

          private:
              
            using ZeroLagXcorr<ValueType>::numDimension;
            using ZeroLagXcorr<ValueType>::equationType;
            using ZeroLagXcorr<ValueType>::gradientKernel;
            using ZeroLagXcorr<ValueType>::gradientDomain;
            using ZeroLagXcorr<ValueType>::useSourceEncode;
            using ZeroLagXcorr<ValueType>::decomposition;
            using ZeroLagXcorr<ValueType>::NT;
             
            
            /* required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrRho;
            using ZeroLagXcorr<ValueType>::xcorrLambda;
            using ZeroLagXcorr<ValueType>::xcorrRhoSuRu;
            using ZeroLagXcorr<ValueType>::xcorrRhoSdRd;
            using ZeroLagXcorr<ValueType>::xcorrRhoSuRd;
            using ZeroLagXcorr<ValueType>::xcorrRhoSdRu;
            using ZeroLagXcorr<ValueType>::xcorrLambdaSuRu;
            using ZeroLagXcorr<ValueType>::xcorrLambdaSdRd;
            using ZeroLagXcorr<ValueType>::xcorrLambdaSuRd;
            using ZeroLagXcorr<ValueType>::xcorrLambdaSdRu;
            /* non required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrMuA;
            using ZeroLagXcorr<ValueType>::xcorrMuB;
            using ZeroLagXcorr<ValueType>::xcorrMuC;
            
            std::string type;
        };
    }
}
