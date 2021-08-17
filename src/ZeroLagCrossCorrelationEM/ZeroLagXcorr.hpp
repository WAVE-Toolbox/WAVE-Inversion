

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <Configuration/Configuration.hpp>
#include <Wavefields/WavefieldsFactory.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "../WorkflowEM/Workflow.hpp"
#include "../Common/Common.hpp"
#include "../Taper/Taper2D.hpp"

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
        class ZeroLagXcorrEM
        {

          public:
            //! Default constructor
            ZeroLagXcorrEM(){};
            //! Default destructor
            ~ZeroLagXcorrEM(){};

            //! \brief Declare ZeroLagXcorr pointer
            typedef std::shared_ptr<ZeroLagXcorrEM<ValueType>> ZeroLagXcorrPtr;

            //! Reset cross correlated wavefields
            virtual void resetXcorr(KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) = 0;

            virtual void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) = 0;
            virtual void updateHessianVectorProduct(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefieldDerivative, Wavefields::Wavefields<ValueType> &adjointWavefield, Wavefields::Wavefields<ValueType> &forwardWavefield2ndOrder, Wavefields::Wavefields<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) = 0;

            virtual int getNumDimension() const = 0;
            virtual std::string getEquationType() const = 0;

            virtual scai::lama::DenseVector<ValueType> const &getXcorrSigmaEM() const;  
            virtual scai::lama::DenseVector<ValueType> const &getXcorrEpsilonEM() const; 
            virtual scai::lama::DenseVector<ValueType> const &getXcorrRSigmaEM() const;   
            virtual scai::lama::DenseVector<ValueType> const &getXcorrREpsilonEM() const;

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) = 0;

            virtual void write(std::string filename, scai::IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) = 0;
            void setDecomposeType(scai::IndexType setDecomposeType);
            void setGradientType(scai::IndexType setGradientType);

          protected:
            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            void writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string type, scai::IndexType t);

            int numDimension;
            std::string equationType; 

            scai::lama::DenseVector<ValueType> xcorrSigmaEM; //!< correlated Wavefields for the sigma gradientEM
            scai::lama::DenseVector<ValueType> xcorrEpsilonEM; //!< correlated Wavefields for the epsilon gradientEM
            scai::lama::DenseVector<ValueType> xcorrRSigmaEM; //!< correlated Wavefields for the rSigma gradientEM
            scai::lama::DenseVector<ValueType> xcorrREpsilonEM;    //!< correlated Wavefields for the rEpsilonEM gradientEM
            scai::lama::DenseVector<ValueType> xcorrSigmaEMstep; //!< correlated Wavefields for the sigma gradientEM
            scai::lama::DenseVector<ValueType> xcorrEpsilonEMstep; //!< correlated Wavefields for the epsilon gradientEM
            scai::lama::DenseVector<ValueType> xcorrSigmaEMSuRu;
            scai::lama::DenseVector<ValueType> xcorrSigmaEMSdRd;
            scai::lama::DenseVector<ValueType> xcorrSigmaEMSuRd;
            scai::lama::DenseVector<ValueType> xcorrSigmaEMSdRu;
            scai::lama::DenseVector<ValueType> xcorrEpsilonEMSuRu;
            scai::lama::DenseVector<ValueType> xcorrEpsilonEMSdRd;
            scai::lama::DenseVector<ValueType> xcorrEpsilonEMSuRd;
            scai::lama::DenseVector<ValueType> xcorrEpsilonEMSdRu;
            scai::lama::DenseVector<ValueType> xcorrRSigmaEMstep; //!< correlated Wavefields for the rSigma gradientEM
            scai::lama::DenseVector<ValueType> xcorrREpsilonEMstep;    //!< correlated Wavefields for the rEpsilonEM gradientEM
            scai::IndexType decomposeType = 0;
            scai::IndexType gradientType = 0;
        };
    }
}
