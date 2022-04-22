#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/common/Walltime.hpp>

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
            virtual void gatherWavefields(Wavefields::Wavefields<ValueType> &wavefields, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint) = 0;
            virtual void sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sourceFC, ValueType DT, scai::IndexType shotNumber, std::vector<scai::lama::SparseVector<ValueType>> taperEncode) = 0;
            virtual void applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

            virtual int getNumDimension() const = 0;
            virtual std::string getEquationType() const = 0;

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot) = 0;

            virtual void write(std::string filename, scai::IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;
            void prepareForInversion(scai::IndexType setGradientKernel, KITGPI::Configuration::Configuration config);

            /* Seismic */
            virtual scai::lama::DenseVector<ValueType> const &getXcorrRho() const = 0;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrLambda() const = 0;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuA() const = 0;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuB() const = 0;
            virtual scai::lama::DenseVector<ValueType> const &getXcorrMuC() const = 0;
            
            /* EM */
            virtual scai::lama::DenseVector<ValueType> const &getXcorrSigma() const = 0;  
            virtual scai::lama::DenseVector<ValueType> const &getXcorrEpsilon() const = 0; 
            virtual scai::lama::DenseVector<ValueType> const &getXcorrREpsilonSigma() const = 0; 

          protected:
            /* Common */
            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            void writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string type, scai::IndexType t);

            typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;
            int numDimension;
            std::string equationType; 
            
            scai::IndexType NT;
            scai::IndexType decomposition = 0;
            scai::IndexType gradientKernel = 0;
            scai::IndexType gradientDomain = 0;
            bool normalizeGradient = false;
            scai::IndexType useSourceEncode = 0;
            scai::IndexType numRelaxationMechanisms = 0; //!< Number of relaxation mechanisms
            std::vector<ValueType> relaxationFrequency; 
            
            /* Seismic */
            scai::lama::DenseVector<ValueType> xcorrMuA; //!< correlated Wavefields for the Mu gradient
            scai::lama::DenseVector<ValueType> xcorrMuB; //!< correlated Wavefields for the Mu gradient
            scai::lama::DenseVector<ValueType> xcorrMuC; //!< correlated Wavefields for the Mu gradient

            scai::lama::DenseVector<ValueType> xcorrRho;    //!< correlated Wavefields for the rho gradient
            scai::lama::DenseVector<ValueType> xcorrLambda; //!< correlated Wavefields for the lambda gradient
            scai::lama::DenseVector<ValueType> xcorrRhoSuRu; 
            scai::lama::DenseVector<ValueType> xcorrRhoSdRd; 
            scai::lama::DenseVector<ValueType> xcorrRhoSuRd;
            scai::lama::DenseVector<ValueType> xcorrRhoSdRu; 
            scai::lama::DenseVector<ValueType> xcorrLambdaSuRu; 
            scai::lama::DenseVector<ValueType> xcorrLambdaSdRd;
            scai::lama::DenseVector<ValueType> xcorrLambdaSuRd;
            scai::lama::DenseVector<ValueType> xcorrLambdaSdRu;
            scai::lama::DenseMatrix<ValueType> VXforward;
            scai::lama::DenseMatrix<ValueType> VXadjoint;
            scai::lama::DenseMatrix<ValueType> VYforward;
            scai::lama::DenseMatrix<ValueType> VYadjoint;
            scai::lama::DenseMatrix<ValueType> VZforward;
            scai::lama::DenseMatrix<ValueType> VZadjoint;
            scai::lama::DenseMatrix<ValueType> Sxxforward;
            scai::lama::DenseMatrix<ValueType> Sxxadjoint;
            scai::lama::DenseMatrix<ValueType> Syyforward;
            scai::lama::DenseMatrix<ValueType> Syyadjoint;
            scai::lama::DenseMatrix<ValueType> Sxyforward;
            scai::lama::DenseMatrix<ValueType> Sxyadjoint;
            scai::lama::DenseMatrix<ValueType> Sxzforward;
            scai::lama::DenseMatrix<ValueType> Sxzadjoint;
            scai::lama::DenseMatrix<ValueType> Syzforward;
            scai::lama::DenseMatrix<ValueType> Syzadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fVXforward;
            scai::lama::DenseMatrix<ComplexValueType> fVXadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fVYforward;
            scai::lama::DenseMatrix<ComplexValueType> fVYadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fVZforward;
            scai::lama::DenseMatrix<ComplexValueType> fVZadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fSxxforward;
            scai::lama::DenseMatrix<ComplexValueType> fSxxadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fSyyforward;
            scai::lama::DenseMatrix<ComplexValueType> fSyyadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fSxyforward;
            scai::lama::DenseMatrix<ComplexValueType> fSxyadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fSxzforward;
            scai::lama::DenseMatrix<ComplexValueType> fSxzadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fSyzforward;
            scai::lama::DenseMatrix<ComplexValueType> fSyzadjoint;
            
            /* EM */
            scai::lama::DenseVector<ValueType> xcorrSigma; //!< correlated Wavefields for the sigma gradient
            scai::lama::DenseVector<ValueType> xcorrEpsilon; //!< correlated Wavefields for the epsilon gradient
            scai::lama::DenseVector<ValueType> xcorrREpsilonSigma; //!< correlated Wavefields for the rEpsilonSigma gradient
            scai::lama::DenseVector<ValueType> xcorrSigmaSuRu;
            scai::lama::DenseVector<ValueType> xcorrSigmaSdRd;
            scai::lama::DenseVector<ValueType> xcorrSigmaSuRd;
            scai::lama::DenseVector<ValueType> xcorrSigmaSdRu;
            scai::lama::DenseVector<ValueType> xcorrEpsilonSuRu;
            scai::lama::DenseVector<ValueType> xcorrEpsilonSdRd;
            scai::lama::DenseVector<ValueType> xcorrEpsilonSuRd;
            scai::lama::DenseVector<ValueType> xcorrEpsilonSdRu;
            scai::lama::DenseMatrix<ValueType> EXforward;
            scai::lama::DenseMatrix<ValueType> EXadjoint;
            scai::lama::DenseMatrix<ValueType> EYforward;
            scai::lama::DenseMatrix<ValueType> EYadjoint;
            scai::lama::DenseMatrix<ValueType> EZforward;
            scai::lama::DenseMatrix<ValueType> EZadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fEXforward;
            scai::lama::DenseMatrix<ComplexValueType> fEXadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fEYforward;
            scai::lama::DenseMatrix<ComplexValueType> fEYadjoint;
            scai::lama::DenseMatrix<ComplexValueType> fEZforward;
            scai::lama::DenseMatrix<ComplexValueType> fEZadjoint;
        };
    }
}
