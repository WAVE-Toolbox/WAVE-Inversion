

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

            //! Reset cross correlated wavefields
            virtual void resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

            virtual void update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

            virtual int getNumDimension() const = 0;
            virtual std::string getEquationType() const = 0;

            scai::lama::DenseVector<ValueType> const &getXcorrRho() const;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrLambda() const;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuA() const;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuB() const;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuC() const;

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

            virtual void write(std::string type, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

          protected:
            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            void writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string type, scai::IndexType t);

            int numDimension;
            std::string equationType; 

            scai::lama::DenseVector<ValueType> xcorrMuA; //!< correlated Wavefields for the Mu gradient
            scai::lama::DenseVector<ValueType> xcorrMuB; //!< correlated Wavefields for the Mu gradient
            scai::lama::DenseVector<ValueType> xcorrMuC; //!< correlated Wavefields for the Mu gradient

            scai::lama::DenseVector<ValueType> xcorrRho;    //!< correlated Wavefields for the rho gradient
            scai::lama::DenseVector<ValueType> xcorrLambda; //!< correlated Wavefields for the lambda gradient
        };
    }
}
