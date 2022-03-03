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
        explicit SourceEstimation() : useOffsetMutes(false), readTaper(false), taperName(""){};
        ~SourceEstimation(){};

        void init(scai::IndexType nt, scai::dmemo::DistributionPtr sourceDistribution, ValueType waterLvl, std::string tprName = "");
        void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr sourceDistribution, Taper::Taper1D<ValueType> &sourceSignalTaper);

        typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;

        void estimateSourceSignal(KITGPI::Acquisition::Receivers<ValueType> const &receivers, KITGPI::Acquisition::Receivers<ValueType> const &receiversTrue, scai::IndexType shotInd, scai::IndexType shotNumber);
        void applyFilter(KITGPI::Acquisition::Sources<ValueType> &sourcesEncode, IndexType shotNumberEncode, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode) const;
        
        void estimateSourceSignalEncode(scai::dmemo::CommunicatorPtr commShot, scai::IndexType shotNumberEncode, KITGPI::Configuration::Configuration const &config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, KITGPI::Acquisition::Receivers<ValueType> const &receivers, KITGPI::Acquisition::Receivers<ValueType> const &receiversTrue);
        
        void calcOffsetMutes(KITGPI::Acquisition::Sources<ValueType> const &sources, KITGPI::Acquisition::Receivers<ValueType> &receivers, ValueType minOffset, ValueType maxOffset, scai::IndexType shotInd, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates);
        void calcRefTraces(Configuration::Configuration const &config, scai::IndexType shotInd, KITGPI::Acquisition::Receivers<ValueType> &receivers, Taper::Taper1D<ValueType> const &sourceSignalTaper);
        void setRefTracesToSource(KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> const &receivers, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, scai::IndexType shotIndTrue, scai::IndexType shotNumberEncode);
        void calcOffsetMutesEncode(scai::dmemo::CommunicatorPtr commShot, scai::IndexType shotNumberEncode, KITGPI::Configuration::Configuration const &config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, KITGPI::Acquisition::Receivers<ValueType> &receiversEncode);
        void calcRefTracesEncode(scai::dmemo::CommunicatorPtr commShot, scai::IndexType shotNumberEncode, KITGPI::Configuration::Configuration const &config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, KITGPI::Acquisition::Receivers<ValueType> &receiversEncode, Taper::Taper1D<ValueType> const &sourceSignalTaper);
        
      private:
        ValueType waterLevel;
        scai::IndexType nFFT; // filter length
        bool isSeismic;
        scai::IndexType useSourceEncode;

        scai::lama::DenseMatrix<ComplexValueType> filter;

        bool useOffsetMutes;
        std::vector<scai::lama::DenseVector<ValueType>> mutes;
        std::vector<scai::lama::DenseVector<ValueType>> offsets;
        std::vector<scai::lama::DenseMatrix<ValueType>> refTraces;
        bool readTaper;
        std::string taperName;

        void matCorr(scai::lama::DenseVector<ComplexValueType> &prod, scai::lama::DenseMatrix<ValueType> const &A, scai::lama::DenseMatrix<ValueType> const &B, scai::IndexType shotInd);
        void addComponents(scai::lama::DenseVector<ComplexValueType> &sum, KITGPI::Acquisition::Receivers<ValueType> const &receiversA, KITGPI::Acquisition::Receivers<ValueType> const &receiversB, scai::IndexType shotInd, scai::IndexType shotNumber);
    };
}
