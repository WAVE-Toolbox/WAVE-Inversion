#pragma once

#include <scai/lama.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>

#include <Configuration/Configuration.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>

#include "../Gradient/GradientFactory.hpp"
#include "Misfit/Misfit.hpp"
#include "Misfit/MisfitFactory.hpp"

/*! \brief The class StepLengthSearch searches for the optimal steplength for model update
 *
 */
template <typename ValueType>
class StepLengthSearch{
    
public:
    
    StepLengthSearch() : step2ok(false), step3ok(false), stepCalcCount(0), steplengthOptimum(0), steplengthParabola(3), misfitParabola(3) {};
    ~StepLengthSearch(){};
    
    void run(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, scai::lama::Scalar steplength_init, scai::lama::Scalar currentMisfit);
    
    void initLogFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, std::string misfitType);
    void appendToLogFile(scai::dmemo::CommunicatorPtr comm, IndexType iteration, std::string logFilename, scai::lama::Scalar misfitSum);
    
    scai::lama::Scalar const &getSteplength();
    scai::lama::Scalar parabolicFit(scai::lama::DenseVector<ValueType> const &steplengthParabola,scai::lama::DenseVector<ValueType> const &misfitParabola);
    
private:
    
    scai::lama::Scalar calcMisfit(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, scai::lama::Scalar steplength);
    
    bool step2ok;
    bool step3ok;
    int stepCalcCount;
    scai::lama::Scalar steplengthOptimum;
    scai::lama::Scalar steplengthMin;
    scai::lama::Scalar steplengthMax;
    scai::lama::DenseVector<ValueType> steplengthParabola;
    scai::lama::DenseVector<ValueType> misfitParabola;
    //scai::lama::Scalar steplengthExtremum;

    typedef typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
    wavefieldPtr wavefields;
    
    std::ofstream logFile;

    
};
