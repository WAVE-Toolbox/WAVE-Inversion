

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

        /*! \brief Class to hold and caclulate the zero lag cross correlated wavefields for 2D elastic gradients 
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Delastic : public ZeroLagXcorrSeismic<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Delastic(){equationType="elastic"; numDimension=2;};

            //! Default destructor
            ~ZeroLagXcorr2Delastic(){};

            explicit ZeroLagXcorr2Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot);

            void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void gatherWavefields(Wavefields::Wavefields<ValueType> &wavefields, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint) override;
            void sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sourceFC, ValueType DT, scai::IndexType shotNumber) override;
            void applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;

            int getNumDimension() const;
            std::string getEquationType() const;

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
            using ZeroLagXcorr<ValueType>::xcorrRho;
            using ZeroLagXcorr<ValueType>::xcorrMuA;
            using ZeroLagXcorr<ValueType>::xcorrMuB;
            using ZeroLagXcorr<ValueType>::xcorrMuC;
            using ZeroLagXcorr<ValueType>::xcorrLambda;
            using ZeroLagXcorr<ValueType>::VXforward;
            using ZeroLagXcorr<ValueType>::VXadjoint;
            using ZeroLagXcorr<ValueType>::VYforward;
            using ZeroLagXcorr<ValueType>::VYadjoint;
            using ZeroLagXcorr<ValueType>::Sxxforward;
            using ZeroLagXcorr<ValueType>::Sxxadjoint;
            using ZeroLagXcorr<ValueType>::Syyforward;
            using ZeroLagXcorr<ValueType>::Syyadjoint;
            using ZeroLagXcorr<ValueType>::Sxyforward;
            using ZeroLagXcorr<ValueType>::Sxyadjoint;
            using ZeroLagXcorr<ValueType>::fVXforward;
            using ZeroLagXcorr<ValueType>::fVXadjoint;
            using ZeroLagXcorr<ValueType>::fVYforward;
            using ZeroLagXcorr<ValueType>::fVYadjoint;
            using ZeroLagXcorr<ValueType>::fSxxforward;
            using ZeroLagXcorr<ValueType>::fSxxadjoint;
            using ZeroLagXcorr<ValueType>::fSyyforward;
            using ZeroLagXcorr<ValueType>::fSyyadjoint;
            using ZeroLagXcorr<ValueType>::fSxyforward;
            using ZeroLagXcorr<ValueType>::fSxyadjoint;

            std::string type;
        };
    }
}
