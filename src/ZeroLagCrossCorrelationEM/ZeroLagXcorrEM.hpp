
#pragma once

#include "../ZeroLagCrossCorrelation/ZeroLagXcorr.hpp"

namespace KITGPI
{

    //! \brief ZeroLagXcorr namespace
    namespace ZeroLagXcorr
    {

        /*! \brief Abstract class to calculate the cross correlation between forward and adjoint wavefield.
         *
         * ZeroLagXcorr implements some methods, which are required by all derived classes.
         * As this class is an abstract class, all methods are protected.
         */
        template <typename ValueType>
        class ZeroLagXcorrEM : public ZeroLagXcorr<ValueType>
        {

          public:
            //! Default constructor
            ZeroLagXcorrEM(){};
            //! Default destructor
            ~ZeroLagXcorrEM(){};

            //! Reset cross correlated wavefields
            virtual void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

            virtual void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

            virtual int getNumDimension() const = 0;
            virtual std::string getEquationType() const = 0;

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

            virtual void write(std::string filename, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;
            
            /* Seismic */
            virtual scai::lama::DenseVector<ValueType> const &getXcorrRho() const override;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrLambda() const override;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuA() const override;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuB() const override;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuC() const override;
            
            /* EM */
            virtual scai::lama::DenseVector<ValueType> const &getXcorrSigmaEM() const override;  
            virtual scai::lama::DenseVector<ValueType> const &getXcorrEpsilonEM() const override; 
            virtual scai::lama::DenseVector<ValueType> const &getXcorrREpsilonSigmaEM() const override; 

          protected:
            /* Common */
            using ZeroLagXcorr<ValueType>::numDimension;
            using ZeroLagXcorr<ValueType>::equationType; 
            
            using ZeroLagXcorr<ValueType>::decomposeWavefieldType;
            using ZeroLagXcorr<ValueType>::gradientType;

            /* Seismic */
            using ZeroLagXcorr<ValueType>::xcorrMuA; //!< correlated Wavefields for the Mu gradient
            using ZeroLagXcorr<ValueType>::xcorrMuB; //!< correlated Wavefields for the Mu gradient
            using ZeroLagXcorr<ValueType>::xcorrMuC; //!< correlated Wavefields for the Mu gradient

            using ZeroLagXcorr<ValueType>::xcorrRho;    //!< correlated Wavefields for the rho gradient
            using ZeroLagXcorr<ValueType>::xcorrLambda; //!< correlated Wavefields for the lambda gradient
            using ZeroLagXcorr<ValueType>::xcorrRhostep; 
            using ZeroLagXcorr<ValueType>::xcorrLambdastep;
            using ZeroLagXcorr<ValueType>::xcorrRhoSuRu; 
            using ZeroLagXcorr<ValueType>::xcorrRhoSdRd; 
            using ZeroLagXcorr<ValueType>::xcorrRhoSuRd;
            using ZeroLagXcorr<ValueType>::xcorrRhoSdRu; 
            using ZeroLagXcorr<ValueType>::xcorrLambdaSuRu; 
            using ZeroLagXcorr<ValueType>::xcorrLambdaSdRd;
            using ZeroLagXcorr<ValueType>::xcorrLambdaSuRd;
            using ZeroLagXcorr<ValueType>::xcorrLambdaSdRu;
            
            /* EM */
            using ZeroLagXcorr<ValueType>::xcorrSigmaEM; //!< correlated Wavefields for the sigma gradient
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEM; //!< correlated Wavefields for the epsilon gradient
            using ZeroLagXcorr<ValueType>::xcorrREpsilonSigmaEM; //!< correlated Wavefields for the rSigma gradient
            using ZeroLagXcorr<ValueType>::xcorrSigmaEMstep; //!< correlated Wavefields for the sigma gradient
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEMstep; //!< correlated Wavefields for the epsilon gradient
            using ZeroLagXcorr<ValueType>::xcorrSigmaEMSuRu;
            using ZeroLagXcorr<ValueType>::xcorrSigmaEMSdRd;
            using ZeroLagXcorr<ValueType>::xcorrSigmaEMSuRd;
            using ZeroLagXcorr<ValueType>::xcorrSigmaEMSdRu;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEMSuRu;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEMSdRd;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEMSuRd;
            using ZeroLagXcorr<ValueType>::xcorrEpsilonEMSdRu;
            using ZeroLagXcorr<ValueType>::xcorrREpsilonSigmaEMstep; //!< correlated Wavefields for the rSigma gradient gradient
        };
    }
}
