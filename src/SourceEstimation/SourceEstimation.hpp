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
        SourceEstimation(scai::IndexType nt);
        SourceEstimation(ValueType waterLvl, scai::IndexType nt);

        void setWaterLevel(ValueType waterLvl);

        ~SourceEstimation(){};

        typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;

        void estimateSourceSignal(KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Acquisition::Sources<ValueType> &sources);

      private:
        ValueType waterLevel;
        scai::IndexType nFFT; // filter length
        scai::lama::DenseVector<ComplexValueType> filter;

        void matCorr(scai::lama::DenseVector<ComplexValueType> &prod, scai::lama::DenseMatrix<ValueType> const &A, scai::lama::DenseMatrix<ValueType> const &B);
        void addComponents(scai::lama::DenseVector<ComplexValueType> &sum, KITGPI::Acquisition::Receivers<ValueType> const &receiversA, KITGPI::Acquisition::Receivers<ValueType> const &receiversB);
        void applyFilter(KITGPI::Acquisition::Sources<ValueType> &sources);
    };
}