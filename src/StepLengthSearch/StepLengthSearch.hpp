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
#include "../Misfit/Misfit.hpp"
#include "../Misfit/MisfitFactory.hpp"

namespace KITGPI
{
    
    /*! \brief Class to do an inexact line search for finding an optimal steplength for the model update
     * 
     * The inexact line search is done by applying a parabolic fit if appropriate (steplength, misfit) pairs can be found. 
     *
     */
    template <typename ValueType>
    class StepLengthSearch{
        
    public:
        
        StepLengthSearch() : step2ok(false), step3ok(false), stepCalcCount(0), steplengthOptimum(0), steplengthParabola(3, 0), misfitParabola(3, 0) {};
        ~StepLengthSearch(){};
        
        
        void run(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, ValueType steplengthInit, scai::lama::DenseVector<ValueType> currentMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow);
        
        void initLogFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, std::string misfitType);
        void appendToLogFile(scai::dmemo::CommunicatorPtr comm, scai::IndexType workflowStage, scai::IndexType iteration, std::string logFilename, ValueType misfitSum);
        
        ValueType const &getSteplength();
        ValueType parabolicFit(scai::lama::DenseVector<ValueType> const &steplengthParabola,scai::lama::DenseVector<ValueType> const &misfitParabola);
        
    private:
        
        ValueType calcMisfit(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, ValueType steplength, KITGPI::Workflow::Workflow<ValueType> const &workflow);
        
        bool step2ok;
        bool step3ok;
        int stepCalcCount;
        ValueType steplengthOptimum;
        ValueType steplengthMin;
        ValueType steplengthMax;
        scai::lama::DenseVector<ValueType> steplengthParabola;
        scai::lama::DenseVector<ValueType> misfitParabola;
        //ValueType steplengthExtremum;

        typedef typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
        wavefieldPtr wavefields;
        
        std::ofstream logFile;

        
    };
}
