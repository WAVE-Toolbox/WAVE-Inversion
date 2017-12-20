#pragma once

#include <scai/lama.hpp>
#include <scai/common/Walltime.hpp>

#include <Configuration/Configuration.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>

template <typename ValueType>
class StepLengthSearch{
    
public:
    
    StepLengthSearch() : step2ok(false), step3ok(false){};
    ~StepLengthSearch(){};
    
    void calc(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, scai::lama::DenseVector<ValueType> const &update, scai::lama::Scalar steplength_init, scai::lama::Scalar currentMisfit);
    
     scai::lama::Scalar steplengthOptimum; // only for test
    
private:
    
    void parabolicFit();
    scai::lama::Scalar calcMisfit(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, scai::lama::DenseVector<ValueType> const &update, scai::lama::Scalar steplength);
    
    scai::lama::DenseVector<ValueType> steplengthParabola;
    scai::lama::DenseVector<ValueType> misfitParabola;
//     scai::lama::Scalar steplengthParabola1;
//     scai::lama::Scalar steplengthParabola2;
//     scai::lama::Scalar steplengthParabola3;
//     scai::lama::Scalar misfitParabola1;
//     scai::lama::Scalar misfitParabola2;
//     scai::lama::Scalar misfitParabola3;
    scai::lama::Scalar steplengthExtremum;
//     scai::lama::Scalar steplengthOptimum;
    scai::lama::Scalar steplengthMin;
    scai::lama::Scalar steplengthMax;
    bool step2ok;
    bool step3ok;
    
    
};
