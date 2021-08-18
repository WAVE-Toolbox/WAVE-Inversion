#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <Configuration/Configuration.hpp>
#include <Wavefields/WavefieldsFactory.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "../Workflow/Workflow.hpp"

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
        class ZeroLagXcorr
        {

          public:
            //! Default constructor
            ZeroLagXcorr(){};
            //! Default destructor
            ~ZeroLagXcorr(){};

            //! \brief Declare ZeroLagXcorr pointer
            typedef std::shared_ptr<ZeroLagXcorr<ValueType>> ZeroLagXcorrPtr;

            /* Common */
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
            void setDecomposeType(scai::IndexType setDecomposeType);
            void setGradientType(scai::IndexType setGradientType);

            /* Seismic */
            virtual scai::lama::DenseVector<ValueType> const &getXcorrRho() const = 0;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrLambda() const = 0;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuA() const = 0;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuB() const = 0;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuC() const = 0;
            
            /* EM */
            virtual scai::lama::DenseVector<ValueType> const &getXcorrSigmaEM() const = 0;  
            virtual scai::lama::DenseVector<ValueType> const &getXcorrEpsilonEM() const = 0; 
            virtual scai::lama::DenseVector<ValueType> const &getXcorrRSigmaEM() const = 0;   
            virtual scai::lama::DenseVector<ValueType> const &getXcorrREpsilonEM() const = 0;

          protected:
            /* Common */
            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            void writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string type, scai::IndexType t);

            int numDimension;
            std::string equationType; 
            
            scai::IndexType decomposeType = 0;
            scai::IndexType gradientType = 0;

            /* Seismic */
            scai::lama::DenseVector<ValueType> xcorrMuA; //!< correlated Wavefields for the Mu gradient
            scai::lama::DenseVector<ValueType> xcorrMuB; //!< correlated Wavefields for the Mu gradient
            scai::lama::DenseVector<ValueType> xcorrMuC; //!< correlated Wavefields for the Mu gradient

            scai::lama::DenseVector<ValueType> xcorrRho;    //!< correlated Wavefields for the rho gradient
            scai::lama::DenseVector<ValueType> xcorrLambda; //!< correlated Wavefields for the lambda gradient
            scai::lama::DenseVector<ValueType> xcorrRhostep; 
            scai::lama::DenseVector<ValueType> xcorrLambdastep;
            scai::lama::DenseVector<ValueType> xcorrRhoSuRu; 
            scai::lama::DenseVector<ValueType> xcorrRhoSdRd; 
            scai::lama::DenseVector<ValueType> xcorrRhoSuRd;
            scai::lama::DenseVector<ValueType> xcorrRhoSdRu; 
            scai::lama::DenseVector<ValueType> xcorrLambdaSuRu; 
            scai::lama::DenseVector<ValueType> xcorrLambdaSdRd;
            scai::lama::DenseVector<ValueType> xcorrLambdaSuRd;
            scai::lama::DenseVector<ValueType> xcorrLambdaSdRu;
            
            /* EM */
            scai::lama::DenseVector<ValueType> xcorrSigmaEM; //!< correlated Wavefields for the sigma gradient
            scai::lama::DenseVector<ValueType> xcorrEpsilonEM; //!< correlated Wavefields for the epsilon gradient
            scai::lama::DenseVector<ValueType> xcorrRSigmaEM; //!< correlated Wavefields for the rSigma gradient
            scai::lama::DenseVector<ValueType> xcorrREpsilonEM;    //!< correlated Wavefields for the rEpsilonEM gradient
            scai::lama::DenseVector<ValueType> xcorrSigmaEMstep; //!< correlated Wavefields for the sigma gradient
            scai::lama::DenseVector<ValueType> xcorrEpsilonEMstep; //!< correlated Wavefields for the epsilon gradient
            scai::lama::DenseVector<ValueType> xcorrSigmaEMSuRu;
            scai::lama::DenseVector<ValueType> xcorrSigmaEMSdRd;
            scai::lama::DenseVector<ValueType> xcorrSigmaEMSuRd;
            scai::lama::DenseVector<ValueType> xcorrSigmaEMSdRu;
            scai::lama::DenseVector<ValueType> xcorrEpsilonEMSuRu;
            scai::lama::DenseVector<ValueType> xcorrEpsilonEMSdRd;
            scai::lama::DenseVector<ValueType> xcorrEpsilonEMSuRd;
            scai::lama::DenseVector<ValueType> xcorrEpsilonEMSdRu;
            scai::lama::DenseVector<ValueType> xcorrRSigmaEMstep; //!< correlated Wavefields for the rSigma gradient
            scai::lama::DenseVector<ValueType> xcorrREpsilonEMstep;    //!< correlated Wavefields for the rEpsilonEM gradient
        };
    }
}
