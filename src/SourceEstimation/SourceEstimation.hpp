#pragma once

#include <scai/lama.hpp>
#include <scai/common/Walltime.hpp>

#include <Configuration/Configuration.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <Acquisition/SeismogramHandler.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>

#include <scai/common/Complex.hpp>
#include <scai/lama/fft.hpp>
#include <Common/Common.hpp>

namespace KITGPI
{
    
    /*! \brief Class to estimate the source time function from any synthetic source
     * 
     * Wiener filter to minimize the misfit between a synthetic and the observed source signal.
     *
     */
    template <typename ValueType>
    class SourceEstimation{
        
    public: 
        explicit SourceEstimation(ValueType waterLvl) : waterLevel(waterLvl) {};
        ~SourceEstimation(){};
        
        typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;
        
        void estimateSourceSignal(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, scai::IndexType tStepEnd);
        
    private:
        ValueType waterLevel;
        scai::lama::CSRSparseMatrix<ComplexValueType> filter;
        
        void fftMult(scai::lama::DenseVector<ComplexValueType> &prod, scai::lama::DenseMatrix<ValueType> const &A, scai::lama::DenseMatrix<ValueType> const &B);
        void addComponents(scai::lama::DenseVector<ComplexValueType> sum, KITGPI::Acquisition::Receivers<ValueType> &receiversA, KITGPI::Acquisition::Receivers<ValueType> &receiversB);
        void applyFilter(KITGPI::Acquisition::Sources<ValueType> &sources);
    };
}