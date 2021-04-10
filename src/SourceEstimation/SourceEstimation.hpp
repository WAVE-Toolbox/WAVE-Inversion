#pragma once

#include <scai/common/Walltime.hpp>
#include <scai/lama.hpp>

#include <Acquisition/Receivers.hpp>
#include <Acquisition/SeismogramHandler.hpp>
#include <Acquisition/Sources.hpp>
#include <Configuration/Configuration.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>

#include <Common/Common.hpp>
#include <scai/common/Complex.hpp>
#include <scai/lama/fft.hpp>

#include "../Common/Common.hpp"
#include "../Taper/Taper1D.hpp"
#include "../Taper/Taper2D.hpp"

namespace KITGPI
{

    /*! \brief Class to estimate the source time function from any synthetic source
     * 
     * Wiener filter to minimize the misfit between a synthetic and the observed source signal.
     *
     */
    template <typename ValueType>
    class SourceEstimation
    {

      public:
        explicit SourceEstimation() : useOffsetMutes(false), mutes(Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE), readTaper(false), taperName(""){};

        void init(scai::IndexType nt, scai::dmemo::DistributionPtr sourceDistribution, ValueType waterLvl, std::string tprName = "");
        void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr sourceDistribution, Taper::Taper1D<ValueType> &sourceSignalTaper);

        ~SourceEstimation(){};

        typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;

        void estimateSourceSignal(KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, IndexType shotInd, IndexType shotNr);
        void applyFilter(KITGPI::Acquisition::Sources<ValueType> &sources, scai::IndexType shotInd) const;
        void calcOffsetMutes(KITGPI::Acquisition::Sources<ValueType> const &sources, KITGPI::Acquisition::Receivers<ValueType> const &receivers, ValueType maxOffsets, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates);

      private:
        ValueType waterLevel;
        scai::IndexType nFFT; // filter length

        scai::lama::DenseMatrix<ComplexValueType> filter;

        bool useOffsetMutes;
        std::vector<scai::lama::DenseVector<ValueType>> mutes;
        bool readTaper;
        std::string taperName;

        void matCorr(scai::lama::DenseVector<ComplexValueType> &prod, scai::lama::DenseMatrix<ValueType> const &A, scai::lama::DenseMatrix<ValueType> const &B, scai::IndexType iComponent);
        void addComponents(scai::lama::DenseVector<ComplexValueType> &sum, KITGPI::Acquisition::Receivers<ValueType> const &receiversA, KITGPI::Acquisition::Receivers<ValueType> const &receiversB, scai::IndexType shotNumber);
    };
}
