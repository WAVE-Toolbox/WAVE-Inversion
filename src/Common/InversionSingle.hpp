#include <scai/lama.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/dmemo/CommunicatorStack.hpp>

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include <Configuration/Configuration.hpp>
#include <Configuration/ValueType.hpp>
#include <CheckParameter/CheckParameter.hpp>
#include <Filter/Filter.hpp>
#include <Common/Common.hpp>
#include <Partitioning/Partitioning.hpp>

#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolver/ForwardSolverFactory.hpp>
#include <ForwardSolver/SourceReceiverImpl/SourceReceiverImpl.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>

#include "../Misfit/AbortCriterion.hpp"
#include "../Misfit/Misfit.hpp"
#include "../Misfit/MisfitFactory.hpp"
#include "../Optimization/OptimizationFactory.hpp"
#include "../Preconditioning/EnergyPreconditioning.hpp"
#include "../SourceEstimation/SourceEstimation.hpp"
#include "../StepLengthSearch/StepLengthSearch.hpp"
#include "../Taper/Taper1D.hpp"
#include "../Taper/Taper2D.hpp"

#include "../Gradient/GradientCalculation.hpp"
#include "../Gradient/GradientFactory.hpp"
#include "../Workflow/Workflow.hpp"

#include <Common/HostPrint.hpp>

using namespace scai;

namespace KITGPI
{    
    /*! \brief Class to implement the inversion for one wave
     * 
     */
    template <typename ValueType>
    class InversionSingle
    {

    public:
        /* Default constructor and destructor */
        InversionSingle(){};
        ~InversionSingle(){};
        
        InversionSingle(KITGPI::Configuration::Configuration config, IndexType inversionType);
        
        void printConfig(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config, IndexType inversionType, IndexType equationInd);
        
        void init(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config, IndexType inversionType, IndexType equationInd, Acquisition::Coordinates<ValueType> &modelCoordinates, Acquisition::Coordinates<ValueType> &modelCoordinatesBig, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr &dist, scai::dmemo::DistributionPtr &distBig, IndexType maxiterations, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &model, KITGPI::Workflow::Workflow<ValueType> &workflow, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivative, typename ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr &derivativesInversion, KITGPI::StepLengthSearch<ValueType> &SLsearch);
        
        void estimateMemory(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config, IndexType inversionType, IndexType equationInd, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> modelCoordinates, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model);
        
        void initStage(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config, IndexType inversionType, IndexType equationInd, scai::hmemo::ContextPtr ctx, KITGPI::Workflow::Workflow<ValueType> &workflow, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr dataMisfit, bool breakLoop, scai::dmemo::DistributionPtr dist, bool &breakLoopEM);
        
        void calcGradient(scai::dmemo::CommunicatorPtr commAll, scai::dmemo::DistributionPtr dist, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &model, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, KITGPI::Workflow::Workflow<ValueType> &workflow, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfit, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivative, typename ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr &derivativesInversion, KITGPI::StepLengthSearch<ValueType> &SLsearch, Taper::Taper2D<ValueType> modelTaper2DJoint, IndexType maxiterations, IndexType &useRTM, bool &breakLoop, scai::hmemo::ContextPtr ctx, IndexType &seedtime, IndexType inversionType, IndexType equationInd, bool &breakLoopEM, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &modelEM, KITGPI::Configuration::Configuration configEM, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesEM, KITGPI::Workflow::Workflow<ValueType> &workflowEM, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfitEM, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivativeEM, typename ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr &derivativesInversionEM);
        
        void updateModel(scai::dmemo::CommunicatorPtr commAll, scai::dmemo::DistributionPtr dist, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &model, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Workflow::Workflow<ValueType> &workflow, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfit, KITGPI::StepLengthSearch<ValueType> &SLsearch, IndexType &useRTM, bool &breakLoop, IndexType inversionType, IndexType equationInd);
        
        void runExtraModelling(scai::dmemo::CommunicatorPtr commAll, scai::dmemo::DistributionPtr dist, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &model, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, KITGPI::Workflow::Workflow<ValueType> &workflow, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfit, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivative, typename ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr &derivativesInversion, KITGPI::StepLengthSearch<ValueType> &SLsearch, Taper::Taper2D<ValueType> modelTaper2DJoint, IndexType maxiterations, IndexType &useRTM, bool &breakLoop, scai::hmemo::ContextPtr ctx, IndexType &seedtime, IndexType inversionType, IndexType equationInd, bool &breakLoopEM, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &modelEM, KITGPI::Configuration::Configuration configEM, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesEM, KITGPI::Workflow::Workflow<ValueType> &workflowEM, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfitEM, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivativeEM, typename ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr &derivativesInversionEM);        

    private:
        
        double start_t, end_t, start_t_shot, end_t_shot; /* For timing */
        
        std::string dimension;
        std::string equationType;
        bool isSeismic;
        IndexType breakLoopType = 0;
        IndexType exchangeStrategy = 0;
        
        Configuration::Configuration configBig;
        bool useStreamConfig;
        IndexType tStepEnd;
        std::string misfitType;
        std::string multiMisfitType;
        std::string gradname;
        std::string logFilename;
        ValueType steplengthInit;
        std::string optimizationType;
        IndexType numRelaxationMechanisms;
        IndexType useSourceEncode;
        IndexType gradientDomain;
        IndexType numShotDomains;
        ValueType memWavefiledsStorage = 0;
        
        Acquisition::Receivers<ValueType> receivers; 
        ValueType NXPerShot;
        IndexType numShotPerSuperShot = 1;
        
        IndexType shotDomain;
        dmemo::DistributionPtr distInversion = nullptr;
        
        typename ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives;
        typename ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver;
        typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr modelPriori;
        typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr modelPerShot;
        
        typedef typename Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
        wavefieldPtr wavefields;
        wavefieldPtr wavefieldsTemp;
        wavefieldPtr wavefieldsInversion;
        std::vector<wavefieldPtr> wavefieldrecord;
        std::vector<wavefieldPtr> wavefieldrecordReflect;
        
        Acquisition::Sources<ValueType> sources;
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings;
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsEncode;
        std::vector<Acquisition::coordinate3D> cutCoordinates;
        std::vector<IndexType> uniqueShotNos;
        std::vector<IndexType> uniqueShotNosEncode;
        IndexType useRandomSource = 0;
        IndexType numshots = 1;
        std::shared_ptr<const dmemo::BlockDistribution> shotDist;
        IndexType maxcount = 1;  
          
        IndexType gradientKernel; 
        IndexType decomposition; 
        IndexType snapType;
        
        Acquisition::Receivers<ValueType> receiversTrue;
        Acquisition::Receivers<ValueType> receiversStart;
        Acquisition::Receivers<ValueType> adjointSources;
        Acquisition::Receivers<ValueType> sourcesReflect;
        Taper::Taper2D<ValueType> seismogramTaper2D;
        Taper::Taper1D<ValueType> seismogramTaper1D;
        
        lama::DenseVector<ValueType> misfitPerIt;
        ValueType weightingStabilizingFunctionalGradient = 1.0;
        ValueType weightingCrossGradient = 1.0;
        
        AbortCriterion<ValueType> abortCriterion;
        SourceEstimation<ValueType> sourceEst;
        Taper::Taper1D<ValueType> sourceSignalTaper;
        Filter::Filter<ValueType> freqFilter;
        std::string transFcnFmly = "butterworth";
        
        typedef typename Gradient::Gradient<ValueType>::GradientPtr gradientPtr;
        gradientPtr gradient;
        gradientPtr gradientPerShot;
        gradientPtr stabilizingFunctionalGradient;
        GradientCalculation<ValueType> gradientCalculation;
        Taper::Taper1D<ValueType> gradientTaper1D;
        Taper::Taper2D<ValueType> wavefieldTaper2D;
        
        Preconditioning::EnergyPreconditioning<ValueType> energyPrecond;
        Preconditioning::EnergyPreconditioning<ValueType> energyPrecondReflect;
        typename Optimization::Optimization<ValueType>::OptimizationPtr gradientOptimization;
        
        std::vector<IndexType> shotHistory;
        std::vector<IndexType> misfitTypeHistory;
    };
}
