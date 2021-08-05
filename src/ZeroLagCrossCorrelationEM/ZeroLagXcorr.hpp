

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <Configuration/Configuration.hpp>
#include <WavefieldsEM/WavefieldsFactory.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "../WorkflowEM/Workflow.hpp"
#include "../Common/Common.hpp"
#include "../Taper/Taper2D.hpp"
#include <Common/Hilbert.hpp>

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

            virtual void update(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow, scai::IndexType gradientType, scai::IndexType decomposeType) = 0;
            virtual void updateHessianVectorProduct(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, Wavefields::WavefieldsEM<ValueType> &forwardWavefield2ndOrder, Wavefields::WavefieldsEM<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow) = 0;

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
            scai::lama::DenseVector<ValueType> xcorrRSigmaEMstep; //!< correlated Wavefields for the rSigma gradientEM
            scai::lama::DenseVector<ValueType> xcorrREpsilonEMstep;    //!< correlated Wavefields for the rEpsilonEM gradientEM
        };
    }
}
