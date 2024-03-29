

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

        /*! \brief Class to hold and caclulate the zero lag cross correlated wavefields for 2D sh gradients
         *
         */
        template <typename ValueType>
        class ZeroLagXcorr2Dsh : public ZeroLagXcorrSeismic<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorr2Dsh(){equationType="sh"; numDimension=2;};

            //! Default destructor
            ~ZeroLagXcorr2Dsh(){};

            explicit ZeroLagXcorr2Dsh(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot);

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
            scai::lama::DenseVector<ValueType> const &getXcorrLambda() const override;         

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
            using ZeroLagXcorr<ValueType>::xcorrMuC;
            using ZeroLagXcorr<ValueType>::VZforward;
            using ZeroLagXcorr<ValueType>::VZadjoint;
            using ZeroLagXcorr<ValueType>::Sxzforward;
            using ZeroLagXcorr<ValueType>::Sxzadjoint;
            using ZeroLagXcorr<ValueType>::Syzforward;
            using ZeroLagXcorr<ValueType>::Syzadjoint;
            using ZeroLagXcorr<ValueType>::fVZforward;
            using ZeroLagXcorr<ValueType>::fVZadjoint;
            using ZeroLagXcorr<ValueType>::fSxzforward;
            using ZeroLagXcorr<ValueType>::fSxzadjoint;
            using ZeroLagXcorr<ValueType>::fSyzforward;
            using ZeroLagXcorr<ValueType>::fSyzadjoint;
            /* non required wavefields */
            using ZeroLagXcorr<ValueType>::xcorrLambda;
            using ZeroLagXcorr<ValueType>::xcorrMuA;
            using ZeroLagXcorr<ValueType>::xcorrMuB;
            
            std::string type;
        };
    }
}
