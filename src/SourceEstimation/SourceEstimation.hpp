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
        SourceEstimation(ValueType waterLvl, scai::IndexType nt) : waterLevel(waterLvl), tStepEnd(nt) {nFFT = Common::calcNextPowTwo<ValueType>(nt-1);}
        ~SourceEstimation(){};
        
        typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;
        
        void estimateSourceSignal(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives);
        
    private:
        ValueType waterLevel;
        scai::IndexType nFFT; // filter length
        scai::IndexType tStepEnd; // seismogram length + 1
        scai::lama::DenseVector<ComplexValueType> filter;
        
        void matCorr(scai::lama::DenseVector<ComplexValueType> &prod, scai::lama::DenseMatrix<ValueType> const &A, scai::lama::DenseMatrix<ValueType> const &B);
        void addComponents(scai::lama::DenseVector<ComplexValueType> &sum, KITGPI::Acquisition::Receivers<ValueType> const &receiversA, KITGPI::Acquisition::Receivers<ValueType> const &receiversB);
        void applyFilter(KITGPI::Acquisition::Sources<ValueType> &sources);
    };
}