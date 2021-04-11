#include <scai/lama.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/dmemo/CommunicatorStack.hpp>

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include <Configuration/Configuration.hpp>
#include <Configuration/ValueType.hpp>
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>

#include <ForwardSolver/ForwardSolver.hpp>

#include <CheckParameter/CheckParameter.hpp>

#include <Filter/Filter.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolver/ForwardSolverFactory.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>

#include <Partitioning/Partitioning.hpp>

#include <Wavefields/WavefieldsFactory.hpp>

#include "Gradient/GradientCalculation.hpp"
#include "Gradient/GradientFactory.hpp"
#include "Misfit/AbortCriterion.hpp"
#include "Misfit/Misfit.hpp"
#include "Misfit/MisfitFactory.hpp"
#include "Optimization/OptimizationFactory.hpp"
#include "Preconditioning/EnergyPreconditioning.hpp"
#include "SourceEstimation/SourceEstimation.hpp"
#include "StepLengthSearch/StepLengthSearch.hpp"
#include "Taper/Taper1D.hpp"
#include "Taper/Taper2D.hpp"
#include "Workflow/Workflow.hpp"

#include <Common/HostPrint.hpp>

using namespace scai;
using namespace KITGPI;

bool verbose; // global variable definition

int main(int argc, char *argv[])
{
    double start_t, end_t, start_t_shot, end_t_shot; /* For timing */
    double globalStart_t, globalEnd_t; /* For timing */
    globalStart_t = common::Walltime::get();
    
    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    /* --------------------------------------- */
    /* Read configuration from file            */
    /* --------------------------------------- */
    Configuration::Configuration config(argv[1]);
    
    Configuration::Configuration configBig;
    bool useStreamConfig;
    try {
        useStreamConfig = config.get<bool>("useStreamConfig");
    } catch(...) {
        useStreamConfig = false;
    }
    
    if (useStreamConfig){
        configBig.readFromFile(config.get<std::string>("streamConfigFilename"),true);
    }

    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);   
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower); 

    verbose = config.get<bool>("verbose");
    
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    std::string misfitType = config.get<std::string>("misfitType");
    std::string fieldSeisName(config.get<std::string>("fieldSeisName"));
    std::string gradname(config.get<std::string>("gradientFilename"));
    std::string logFilename = config.get<std::string>("logFilename");
    ValueType steplengthInit = config.get<ValueType>("steplengthInit");
    IndexType maxiterations = config.get<IndexType>("maxIterations");
    IndexType maxOutShotIteration =config.get<IndexType>("maxIterations");
    IndexType maxcount = config.get<IndexType>("maxiterations");
    std::string optimizationType = config.get<std::string>("optimizationType");
    
    /* inter node communicator */
    dmemo::CommunicatorPtr commAll = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(commAll->getNodeRank());

    
    HOST_PRINT(commAll, "\n WAVE-Inversion " << dimension << " " << equationType << " - LAMA Version\n");
    HOST_PRINT(commAll, "","  - Running on " << commAll->getSize() << " mpi processes -\n\n");
    
    if (commAll->getRank() == MASTERGPI) {
        config.print();
    }

    std::string settingsFilename; // filename for processor specific settings
    if (common::Settings::getEnvironment(settingsFilename, "SCAI_SETTINGS")) {
        // each processor reads line of settings file that matches its node name and node rank
        common::Settings::readSettingsFile(settingsFilename.c_str(), commAll->getNodeName(), commAll->getNodeRank());
    }
    
    
    /* --------------------------------------- */
    /* coordinate mapping (3D<->1D)            */
    /* --------------------------------------- */

    Acquisition::Coordinates<ValueType> modelCoordinates(config);
    Acquisition::Coordinates<ValueType> modelCoordinatesInversion(config, config.get<IndexType>("DHInversion"));

    if (config.get<bool>("useVariableGrid")) {
        CheckParameter::checkVariableGrid(config, commAll, modelCoordinates);
        for (int layer=0;layer<modelCoordinates.getNumLayers();layer++){
            HOST_PRINT(commAll, "\n Number of gridpoints in layer: " << layer << " = " << modelCoordinates.getNGridpoints(layer));
        }
        auto numGridpointsRegular=config.get<IndexType>("NX")*config.get<IndexType>("NY")*config.get<IndexType>("NZ");
        HOST_PRINT(commAll, "\n Number of gripoints total: " << modelCoordinates.getNGridpoints());
        HOST_PRINT(commAll, "\n Percentage of gridpoints of the underlying regular grid given by NX*NY*NZ: "  << (float) modelCoordinates.getNGridpoints()/numGridpointsRegular * 100 << "% \n\n");
    }

    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */

    
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    IndexType shotDomain = Partitioning::getShotDomain(config, commAll); // will contain the domain to which this processor belongs

    // Build pershots of processors for the shots

    dmemo::CommunicatorPtr commShot = commAll->split(shotDomain);

    dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());
    SCAI_DMEMO_TASK(commShot)

    /* --------------------------------------- */
    /* Distribution                            */
    /* --------------------------------------- */
    
    // distribute the grid onto available processors
    dmemo::DistributionPtr dist = nullptr;
    dmemo::DistributionPtr distInversion = nullptr;
    if ((config.get<IndexType>("partitioning") == 0) || (config.get<IndexType>("partitioning") == 2)) {
        //Block distribution = starting distribution for graph partitioner
        dist = std::make_shared<dmemo::BlockDistribution>(modelCoordinates.getNGridpoints(), commShot);
    } else if (config.get<IndexType>("partitioning") == 1) {
        SCAI_ASSERT(!config.get<bool>("useVariableGrid"), "Grid distribution is not available for the variable grid");
        dist = Partitioning::gridPartition<ValueType>(config, commShot);
        distInversion = KITGPI::Partitioning::gridPartitionInversion<ValueType>(config, commShot);
    } else {
        COMMON_THROWEXCEPTION("unknown partitioning method");
    }
  
    
    if ((config.get<bool>("coordinateWrite")) && (shotDomain==0)) {
        modelCoordinates.writeCoordinates(dist, ctx, config.get<std::string>("coordinateFilename"),config.get<IndexType>("FileFormat"));
    }
    
    /* --------------------------------------- */
    /* Factories                               */
    /* --------------------------------------- */

    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr modelBig(Modelparameter::Factory<ValueType>::Create(equationType));
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(Wavefields::Factory<ValueType>::Create(dimension, equationType));
    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));

    /* --------------------------------------- */
    /* Memory estimation                       */
    /* --------------------------------------- */
    HOST_PRINT(commAll, "============== Memory Estimation: ===============\n\n")

    ValueType memDerivatives = derivatives->estimateMemory(config, dist, modelCoordinates);
    ValueType memWavefileds = wavefields->estimateMemory(dist);
    ValueType memWavefiledsStorage = memWavefileds* (int) (tStepEnd/config.get<IndexType>("DTInversion"));
    ValueType memModel = model->estimateMemory(dist);
    ValueType memSolver = solver->estimateMemory(config, dist, modelCoordinates);
    ValueType memTotal = memDerivatives + memWavefileds + memModel + memSolver + memWavefiledsStorage;

    HOST_PRINT(commAll, " -  Derivative Matrices \t" << memDerivatives << " MB\n");
    HOST_PRINT(commAll, " -  Wavefield vectors \t\t" << memWavefileds << " MB\n");
    HOST_PRINT(commAll, " -  Forward wavefield storage \t" << memWavefiledsStorage << " MB\n");
    HOST_PRINT(commAll, " -  Model Vectors \t\t" << memModel << " MB\n");
    HOST_PRINT(commAll, " -  Boundary Condition Vectors \t" << memSolver << " MB\n");
    HOST_PRINT(commAll, "\n Memory Usage (total / per partition): \n " << memTotal << " / " << memTotal / dist->getNumPartitions() << " MB ");
    IndexType numShotDomains = config.get<IndexType>("NumShotDomains"); // total number of shot domains
    if (numShotDomains > 1)
        HOST_PRINT(commAll, "\n Total Memory Usage (" << numShotDomains << " shot Domains ): \n " << memTotal * numShotDomains << " MB  ");

    HOST_PRINT(commAll, "\n\n========================================================================\n\n")
    
    /* --------------------------------------- */
    /* Call partitioner */
    /* --------------------------------------- */
    if (config.get<IndexType>("partitioning") == 2) {
        start_t = common::Walltime::get();
        dist = Partitioning::graphPartition(config, ctx, commShot, dist, *derivatives,modelCoordinates);
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "", "Finished graph partitioning in " << end_t - start_t << " sec.\n\n");
    }
    
    
    /* --------------------------------------- */
    /* Calculate derivative matrices           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    derivatives->init(dist, ctx, config, modelCoordinates, commShot);
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing matrices in " << end_t - start_t << " sec.\n\n");
    
    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings;
    Acquisition::readAllSettings<ValueType>(sourceSettings, config.get<std::string>("SourceFilename") + ".txt");
    // build Settings for SU?
    //settings = su.getSourceSettings(shotNumber); // currently not working, expecting a sourceSettings struct and not a vector of sourceSettings structs
    //         su.buildAcqMatrixSource(config.get<std::string>("SourceSignalFilename"), modelCoordinates.getDH());
    //         allSettings = su.getSourceSettingsVec();

    Acquisition::Sources<ValueType> sources;

    std::vector<scai::IndexType> uniqueShotNos;
    calcuniqueShotNo(uniqueShotNos, sourceSettings);
    IndexType numshots = uniqueShotNos.size();

    CheckParameter::checkSources<ValueType>(sourceSettings, modelCoordinates, commShot);
    dmemo::BlockDistribution shotDist(numshots, commInterShot);
    
    Acquisition::Receivers<ValueType> receivers;
    if (!config.get<bool>("useReceiversPerShot"))
        receivers.init(config, modelCoordinates, ctx, dist);

    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    // load starting model
    model->init(config, ctx, dist, modelCoordinates);
    //model->init(config, ctx, dist);
    
    /* --------------------------------------- */
    /* Modelparameter for bigmodel             */
    /* --------------------------------------- */
    
    Acquisition::Coordinates<ValueType> modelCoordinatesBig;
    if (useStreamConfig) {
        modelCoordinatesBig.init(configBig);
        dmemo::DistributionPtr distBig = nullptr;
        distBig = std::make_shared<dmemo::BlockDistribution>(modelCoordinatesBig.getNGridpoints(), commShot);
        modelBig->init(configBig, ctx, distBig, modelCoordinatesBig);
    }
    
    std::vector<Acquisition::coordinate3D> cutCoordinates;
    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;
    sources.getAcquisitionSettings(configBig, sourceSettingsBig);
    Acquisition::getCutCoord(cutCoordinates, sourceSettingsBig);
            
    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    wavefields->init(ctx, dist);

    //Temporary Wavefield for the derivative of the forward wavefields
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldsInversion = Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefieldsInversion->init(ctx, distInversion);

    /* --------------------------------------- */
    /* Wavefield record                        */
    /* --------------------------------------- */
    typedef typename Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
    std::vector<wavefieldPtr> wavefieldrecord;

    IndexType dtinversion = config.get<IndexType>("DTInversion");
    for (IndexType i = 0; i < tStepEnd; i++) {
        if (i % dtinversion == 0) {
            wavefieldrecord.push_back(wavefieldsInversion);
        }
    }

    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */
    solver->initForwardSolver(config, *derivatives, *wavefields, *model, modelCoordinates, ctx, config.get<ValueType>("DT"));

    /* --------------------------------------- */
    /* True data                               */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> receiversTrue;
    if (!config.get<bool>("useReceiversPerShot")) {
        receiversTrue.init(config, modelCoordinates, ctx, dist);
    }
    Taper::Taper2D<ValueType> seismogramTaper2D;
    Taper::Taper1D<ValueType> seismogramTaper1D;
    seismogramTaper1D.init(std::make_shared<dmemo::NoDistribution>(tStepEnd), ctx, 1);

    /* --------------------------------------- */
    /* Misfit                                  */
    /* --------------------------------------- */
    Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(Misfit::Factory<ValueType>::Create(misfitType));
    lama::DenseVector<ValueType> misfitPerIt(numshots, 0, ctx);

    /* --------------------------------------- */
    /* Adjoint sources                         */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> adjointSources;
    if (!config.get<bool>("useReceiversPerShot"))
        adjointSources.init(config, modelCoordinates, ctx, dist);

    /* --------------------------------------- */
    /* Workflow                                */
    /* --------------------------------------- */
    Workflow::Workflow<ValueType> workflow(config);

    /* --------------------------------------- */
    /* Abort criterion                         */
    /* --------------------------------------- */
    AbortCriterion<ValueType> abortCriterion;

    /* --------------------------------------- */
    /* Step length search                      */
    /* --------------------------------------- */
    StepLengthSearch<ValueType> SLsearch;
    SLsearch.initLogFile(commAll, logFilename, misfitType);

    /* --------------------------------------- */
    /* Source estimation                       */
    /* --------------------------------------- */
    SourceEstimation<ValueType> sourceEst;
    Taper::Taper1D<ValueType> sourceSignalTaper;
    // calculate source dist
    scai::lama::DenseVector<IndexType> sourcecoords = getsourcecoordinates(sourceSettings, modelCoordinates);
    scai::dmemo::DistributionPtr dist_sources = Acquisition::calcDistribution(sourcecoords, dist);
    if (config.get<bool>("useSourceSignalInversion"))
        sourceEst.init(config, ctx, dist_sources, sourceSignalTaper);
    
    /* --------------------------------------- */
    /* Frequency filter                        */
    /* --------------------------------------- */
    Filter::Filter<ValueType> freqFilter;
    std::string transFcnFmly = "butterworth";
    if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0)
        freqFilter.init(config.get<ValueType>("DT"), tStepEnd);

    /* --------------------------------------- */
    /* Gradients                               */
    /* --------------------------------------- */
    Gradient::Gradient<ValueType>::GradientPtr gradient(Gradient::Factory<ValueType>::Create(equationType));
    gradient->init(ctx, dist);

    Gradient::Gradient<ValueType>::GradientPtr gradientPerShot(Gradient::Factory<ValueType>::Create(equationType));
    gradientPerShot->init(ctx, dist);
    gradientPerShot->setNormalizeGradient(config.get<bool>("normalizeGradient"));
    

//    dmemo::DistributionPtr distBig = nullptr;
//    distBig = std::make_shared<dmemo::BlockDistribution>(modelCoordinatesBig.getNGridpoints(), commShot);
//    Gradient::Gradient<ValueType>::GradientPtr gradientBig(Gradient::Factory<ValueType>::Create(equationType));
//    gradientBig->init(ctx, distBig);


    /* --------------------------------------- */
    /* Gradient calculation                    */
    /* --------------------------------------- */
    GradientCalculation<ValueType> gradientCalculation;

    /* --------------------------------------- */
    /* Gradient taper                          */
    /* --------------------------------------- */
    Preconditioning::SourceReceiverTaper<ValueType> ReceiverTaper;
    if (!config.get<bool>("useReceiversPerShot"))
        ReceiverTaper.init(dist, ctx, receivers, config, modelCoordinates, config.get<IndexType>("receiverTaperRadius"));
    Taper::Taper1D<ValueType> gradientTaper1D;
    if (config.get<bool>("useGradientTaper")) {
        gradientTaper1D.init(dist, ctx, 1);
        gradientTaper1D.read(config.get<std::string>("gradientTaperName"), config.get<IndexType>("FileFormat"));
    }
    bool isSeismic = true;
    Taper::Taper2D<ValueType> gradientTaper2DPerShot;
    gradientTaper2DPerShot.initWavefieldTransform(config, distInversion, dist, ctx, isSeismic); 
    gradientTaper2DPerShot.calcInversionAverageMatrix(modelCoordinates, modelCoordinatesInversion);  

    /* --------------------------------------- */
    /* Gradient preconditioning                */
    /* --------------------------------------- */
    Preconditioning::EnergyPreconditioning<ValueType> energyPrecond;
    if (config.get<bool>("useEnergyPreconditioning") == 1) {
        energyPrecond.init(dist, config);
    }

    /* --------------------------------------- */
    /* Gradient optimization                   */
    /* --------------------------------------- */
    Optimization::Optimization<ValueType>::OptimizationPtr gradientOptimization(Optimization::Factory<ValueType>::Create(optimizationType));
    gradientOptimization->init(dist);

    /* --------------------------------------- */
    /*       Loop over workflow stages         */
    /* --------------------------------------- */

    for (workflow.workflowStage = 0; workflow.workflowStage < workflow.maxStage; workflow.workflowStage++) {

        workflow.printParameters(commAll);

        gradientCalculation.allocate(config, dist, ctx, workflow);

        if (workflow.getLowerCornerFreq() != 0.0 && workflow.getUpperCornerFreq() != 0.0)
            freqFilter.calc(transFcnFmly, "bp", workflow.getFilterOrder(), workflow.getLowerCornerFreq(), workflow.getUpperCornerFreq());
        else if (workflow.getLowerCornerFreq() != 0.0 && workflow.getUpperCornerFreq() == 0.0)
            freqFilter.calc(transFcnFmly, "lp", workflow.getFilterOrder(), workflow.getLowerCornerFreq());
        else if (workflow.getLowerCornerFreq() == 0.0 && workflow.getUpperCornerFreq() != 0.0)
            freqFilter.calc(transFcnFmly, "hp", workflow.getFilterOrder(), workflow.getUpperCornerFreq());
        
        seismogramTaper1D.calcTimeDampingTaper(workflow.getTimeDampingFactor(), config.get<ValueType>("DT"));  
        
        /* --------------------------------------- */
        /*        Loop over pershot shots           */
        /* --------------------------------------- */
        
        std::vector<scai::IndexType> filterHistoryCount(numshots, 0);
        IndexType outShotInd = 0;
        IndexType cutCoordSize = 1;
        if (useStreamConfig) {
            cutCoordSize = cutCoordinates.size();
        }

        while (outShotInd++ < maxOutShotIteration) {
            
            std::srand((int)time(0));
            IndexType cutCoordInd = std::rand() % cutCoordSize;
            while (filterHistoryCount[cutCoordInd] >= maxcount)
                cutCoordInd = std::rand() % cutCoordSize;
            
            if (useStreamConfig==0) {
                outShotInd = maxOutShotIteration;
                cutCoordInd = 0;
            } else {
                modelBig->getModelPerShot(*model, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(cutCoordInd));
            }

            /* --------------------------------------- */
            /*        Loop over iterations             */
            /* --------------------------------------- */
            
            if (useStreamConfig) {
                maxiterations = 1;
            }
            
            for (workflow.iteration = 0; workflow.iteration < maxiterations; workflow.iteration++) {

                HOST_PRINT(commAll, "\n=================================================");
                HOST_PRINT(commAll, "\n============ Workflow stage " << workflow.workflowStage + 1 << " of " << workflow.maxStage << " ==============");
                HOST_PRINT(commAll, "\n============     PerShot " << cutCoordInd + 1  << " of " << cutCoordSize << "     ==============");
                HOST_PRINT(commAll, "\n===========      PerShot Shot " << outShotInd << "      =============");
                HOST_PRINT(commAll, "\n============      Iteration " << workflow.iteration + 1 << "      ==============");
                HOST_PRINT(commAll, "\n=================================================\n\n");
                start_t = common::Walltime::get();
                /* Update model for fd simulation (averaging, inverse Density ...) */
                model->prepareForModelling(modelCoordinates, ctx, dist, commShot);
                
                
                if ((workflow.iteration == 0)&&(commInterShot->getRank() == 0)&&(useStreamConfig == 0)) {
                    /* only shot Domain 0 writes output */
                    model->write((config.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflow.workflowStage+1) + ".It_" + std::to_string(workflow.iteration)), config.get<IndexType>("FileFormat"));
                }
                
                if ((workflow.iteration == 0)&&(useStreamConfig)) {
                    modelBig->write((config.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".pershot_" + std::to_string(cutCoordInd) + ".It_" + std::to_string(workflow.iteration)), config.get<IndexType>("FileFormat"));
                }
                
                solver->prepareForModelling(*model, config.get<ValueType>("DT"));

                /* --------------------------------------- */
                /*        Loop over shots                  */
                /* --------------------------------------- */

                gradient->resetGradient(); // reset gradient because gradient is a sum of all gradientsPerShot gradients+=gradientPerShot
                misfitPerIt = 0;

                IndexType localShotInd = 0;
                IndexType firstShot = shotDist.lb();
                IndexType lastShot = shotDist.ub();
                
                if (useStreamConfig) {
                    lastShot = firstShot + 1;
                }
                
                for (IndexType shotInd = firstShot + cutCoordInd; shotInd < lastShot + cutCoordInd; shotInd++) {
                    IndexType shotNumber = uniqueShotNos[shotInd];
                    localShotInd++;
                    HOST_PRINT(commShot, "Shot number " << shotNumber << ", local shot " << localShotInd << " of " << shotDist.getLocalSize() << ": started\n");

                    if (config.get<bool>("useReceiversPerShot")) {
                        receivers.init(config, modelCoordinates, ctx, dist, shotNumber);
                        receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber);
                        adjointSources.init(config, modelCoordinates, ctx, dist, shotNumber);

                        //CheckParameter::checkReceivers<ValueType>(config, receivers, commShot);

                        ReceiverTaper.init(dist, ctx, receivers, config, modelCoordinates, config.get<IndexType>("receiverTaperRadius"));
                    }

                    /* Read field data (or pseudo-observed data, respectively) */
                    receiversTrue.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), fieldSeisName + ".shot_" + std::to_string(shotNumber), 1);
                                        
                    if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0){
                        receiversTrue.getSeismogramHandler().filter(freqFilter);
                    }
                    
                    if (config.get<IndexType>("useSeismogramTaper") == 2) {
                        seismogramTaper2D.init(receiversTrue.getSeismogramHandler());
                        seismogramTaper2D.read(config.get<std::string>("seismogramTaperName") + ".shot_" + std::to_string(shotNumber) + ".mtx");
                        seismogramTaper2D.apply(receiversTrue.getSeismogramHandler()); 
                    }
                    seismogramTaper1D.apply(receiversTrue.getSeismogramHandler());
                    
                    if (workflow.iteration == 0){
                        receiversTrue.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), fieldSeisName + ".stage_" + std::to_string(workflow.workflowStage+1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
                    }

                    /* Reset approximated Hessian per shot */
                    if (config.get<bool>("useEnergyPreconditioning") == 1)
                        energyPrecond.resetApproxHessian();

                    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
                    Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
                    sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);

                    CheckParameter::checkNumericalArtefeactsAndInstabilities<ValueType>(config, sourceSettingsShot, *model,modelCoordinates,shotNumber);
                    
                    if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0)
                        sources.getSeismogramHandler().filter(freqFilter);
                    
                    /* Source time function inversion */
                    if (config.get<bool>("useSourceSignalInversion")){
                        if (workflow.iteration == 0) {
                            HOST_PRINT(commShot, "Shot number " << shotNumber << ", local shot " << localShotInd << " of " << shotDist.getLocalSize() << " : Source Time Function Inversion\n");

                            wavefields->resetWavefields();

                            for (scai::IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                                solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);
                            }

                            solver->resetCPML();
                            
                            /* Normalize observed and synthetic data */
                            if (config.get<bool>("NormalizeTraces")){
                                receivers.getSeismogramHandler().normalize();
                                receiversTrue.getSeismogramHandler().normalize();
                            }

                            if (config.get<bool>("maxOffsetSrcEst") == 1)
                                sourceEst.calcOffsetMutes(sources, receivers, config.get<ValueType>("maxOffsetSrcEst"),modelCoordinates);
                            
                            sourceEst.estimateSourceSignal(receivers, receiversTrue, shotInd, shotNumber);

                            sourceEst.applyFilter(sources, shotInd);
                            if (config.get<bool>("useSourceSignalTaper"))
                                sourceSignalTaper.apply(sources.getSeismogramHandler());
                            if (config.get<bool>("writeInvertedSource") == 1)
                                sources.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("sourceSeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
                        } else {
                            sourceEst.applyFilter(sources, shotInd);
                            if (config.get<bool>("useSourceSignalTaper"))
                                sourceSignalTaper.apply(sources.getSeismogramHandler());
                        }
                    }
                    
                    /* --------------------------------------- */
                    /*        Forward modelling                */
                    /* --------------------------------------- */

                    HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << numshots << ": Start time stepping with " << tStepEnd << " time steps\n");

                    wavefields->resetWavefields();

                    start_t_shot = common::Walltime::get();
                    //IndexType WavefieldRecordIndex=0;
                    for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {

                        solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);
                        if (tStep % dtinversion == 0) {
                            // save wavefields in std::vector
                            *wavefieldrecord[floor(tStep / dtinversion + 0.5)] =  gradientTaper2DPerShot.applyWavefieldAverage(wavefields);
                        }
                        if (config.get<bool>("useEnergyPreconditioning") == 1) {
                            energyPrecond.intSquaredWavefields(*wavefields, config.get<ValueType>("DT"));
                        }
                    }

                    // check wavefield and seismogram for NaNs or infinite values
                    if ((commShot->any(!wavefields->isFinite(dist)) || commShot->any(!receivers.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write model and break
                    model->write("model_crash", config.get<IndexType>("FileFormat"));
                    COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield, output model as model_crash.FILE_EXTENSION!");
                    }
                    
//                    /* Normalize observed and synthetic data */
//                    if (config.get<bool>("NormalizeTraces")){
//                        receivers.getSeismogramHandler().normalize();
//                        receiversTrue.getSeismogramHandler().normalize();
//                    }

                    if (config.get<IndexType>("useSeismogramTaper") == 2) {                                                   
                        seismogramTaper2D.apply(receivers.getSeismogramHandler()); 
                    }
                    seismogramTaper1D.apply(receivers.getSeismogramHandler());
                    
                    receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration) + ".shot_" + std::to_string(shotNumber), modelCoordinates);

                    HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << numshots << ": Calculate misfit and adjoint sources\n");

                    /* Normalize observed and synthetic data */
                    if (config.get<bool>("NormalizeTraces")){
                        receivers.getSeismogramHandler().normalize();
                        receiversTrue.getSeismogramHandler().normalize();
                    }

                    /* Calculate misfit of one shot */
                    misfitPerIt.setValue(shotInd, dataMisfit->calc(receivers, receiversTrue));

                    /* Calculate adjoint sources */
                    dataMisfit->calcAdjointSources(adjointSources, receivers, receiversTrue);

                    /* Calculate gradient */

                    HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << numshots << ": Start Backward\n");
                    gradientCalculation.run(commAll,*solver, *derivatives, receivers, sources, adjointSources, *model, *gradientPerShot, wavefieldrecord, config, modelCoordinates, shotNumber, workflow, gradientTaper2DPerShot);

                    /* Apply energy preconditioning per shot */
                    if (config.get<bool>("useEnergyPreconditioning") == 1) {
                        energyPrecond.apply(*gradientPerShot, shotNumber, config.get<IndexType>("FileFormat"));
                    }
                    
                    if (config.get<bool>("useReceiversPerShot"))
                        ReceiverTaper.apply(*gradientPerShot);

                    gradientPerShot->normalize();
                    *gradient += *gradientPerShot;
//
//                    if (useStreamConfig) {
//                        IndexType smoothRange = config.get<IndexType>("smoothRange");
//                        IndexType NX = config.get<IndexType>("NX");
//                        IndexType NY = config.get<IndexType>("NY");
//                        IndexType NXBig = configBig.get<IndexType>("NX");
//                        IndexType NYBig = configBig.get<IndexType>("NY");
//                        IndexType boundaryWidth = config.get<IndexType>("BoundaryWidth");
//
//                        gradientBig->setGradientPerShot(*gradientPerShot,modelCoordinates,modelCoordinatesBig,cutCoordinates,cutCoordInd,smoothRange,NX,NY,NXBig,NYBig,boundaryWidth);
//                    }

                    solver->resetCPML();

                    end_t_shot = common::Walltime::get();
                    HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << numshots << ": Finished in " << end_t_shot - start_t_shot << " sec.\n");

                } //end of loop over shots

                gradient->sumShotDomain(commInterShot);
                
                commInterShot->sumArray(misfitPerIt.getLocalValues());

                HOST_PRINT(commAll, "\n======== Finished loop over shots =========");
                HOST_PRINT(commAll, "\n===========================================\n");

                /* Apply receiver Taper (if ReceiverTaperRadius=0 gradient will be multplied by 1) */
                if (!config.get<bool>("useReceiversPerShot"))
                    ReceiverTaper.apply(*gradient);

                gradientOptimization->apply(*gradient, workflow, *model, config);

                if (config.get<IndexType>("FreeSurface") == 2) {
                    lama::DenseVector<ValueType> mask;
                    mask = model->getVelocityP();
                    mask.unaryOp(mask, common::UnaryOp::SIGN);
                    *gradient *= mask;
                }

                if (config.get<bool>("useGradientTaper"))
                    gradientTaper1D.apply(*gradient);

                /* Output of gradient */
                /* only shot Domain 0 writes output */
                if (config.get<IndexType>("WriteGradient") && commInterShot->getRank() == 0) {
                    gradient->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".pershot_" + std::to_string(cutCoordInd) + ".It_" + std::to_string(workflow.iteration + 1), config.get<IndexType>("FileFormat"), workflow);
//                    gradientBig->write(gradnameBig + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".pershot_" + std::to_string(cutCoordInd) + ".It_" + std::to_string(workflow.iteration + 1), config.get<IndexType>("FileFormat"), workflow);
                }
                dataMisfit->addToStorage(misfitPerIt);
                misfitPerIt = 0;

                SLsearch.appendToLogFile(commAll, workflow.workflowStage + 1, workflow.iteration, logFilename, dataMisfit->getMisfitSum(workflow.iteration));

                /* Check abort criteria */
                HOST_PRINT(commAll, "\nMisfit after stage " << workflow.workflowStage + 1 << ", iteration " << workflow.iteration << ": " << dataMisfit->getMisfitSum(workflow.iteration) << "\n");

                bool breakLoop = abortCriterion.check(commAll, *dataMisfit, config, steplengthInit, workflow);
                if (breakLoop == true) {
                    break;
                }

                HOST_PRINT(commAll, "\n===========================================");
                HOST_PRINT(commAll, "\n======== Start step length search =========\n");

                SLsearch.run(commAll, *solver, *derivatives, receivers, sourceSettings, receiversTrue, *model, dist, config, modelCoordinates, *gradient, steplengthInit, dataMisfit->getMisfitIt(workflow.iteration), workflow, freqFilter, sourceEst, sourceSignalTaper, cutCoordInd);

                HOST_PRINT(commAll, "=========== Update Model ============\n\n");
                /* Apply model update */
                *gradient *= SLsearch.getSteplength();
                *model -= *gradient;
                
//                if (useStreamConfig) {
//                    *modelBig -= *gradientBig;
//                }

                if (config.get<bool>("useModelThresholds"))
                    model->applyThresholds(config);

                if ((commInterShot->getRank() == 0)&&(useStreamConfig==0)) {
                    /* only shot Domain 0 writes output */
                        model->write((config.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1)), config.get<IndexType>("FileFormat"));
                }
                
//                 if (useStreamConfig) {
//                     IndexType smoothRange = config.get<IndexType>("smoothRange");
//                     IndexType NXBig = configBig.get<IndexType>("NX");
//                     IndexType NYBig = configBig.get<IndexType>("NY");
//                     IndexType boundaryWidth = config.get<IndexType>("BoundaryWidth");
// 
// //                     modelBig->setModelPerShot(*model,modelCoordinates,modelCoordinatesBig,cutCoordinates,cutCoordInd,smoothRange,NX,NY,NXBig,NYBig,boundaryWidth);
//                     modelBig->write((config.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".pershot_" + std::to_string(cutCoordInd) + ".It_" + std::to_string(workflow.iteration + 1)), config.get<IndexType>("FileFormat"));
//                 }
                
                steplengthInit *= 0.98;

                end_t = common::Walltime::get();
                HOST_PRINT(commAll, "\nFinished iteration " << workflow.iteration + 1 << " in " << end_t - start_t << " sec.\n\n");

                /* -------------------------------------------------------------------- */
                /* One extra forward modelling to ensure complete and consistent output */
                /* -------------------------------------------------------------------- */
                if (workflow.iteration == maxiterations - 1) {

                    HOST_PRINT(commAll, "================ Maximum number of iterations reached =================\n");
                    HOST_PRINT(commAll, "= Do one more forward modelling to calculate misfit and save seismograms =\n\n");

                    /* Update model for fd simulation (averaging, inverse Density ...) */
                    model->prepareForModelling(modelCoordinates, ctx, dist, commShot);
                    solver->prepareForModelling(*model, config.get<ValueType>("DT"));
                    
                    IndexType firstShot = shotDist.lb();
                    IndexType lastShot = shotDist.ub();
                    
                    if (useStreamConfig) {
                        lastShot = firstShot + 1;
                    }
                    
                    for (IndexType shotInd = firstShot + cutCoordInd; shotInd < lastShot + cutCoordInd; shotInd++) {
                        IndexType shotNumber = uniqueShotNos[shotInd];

                        /* Read field data (or pseudo-observed data, respectively) */
                        if (config.get<bool>("useReceiversPerShot")) {
                            receivers.init(config, modelCoordinates, ctx, dist, shotNumber);
                            receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber);
                        }
                        receiversTrue.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), fieldSeisName + ".shot_" + std::to_string(shotNumber), 1);

                        HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << numshots << ": Additional forward run with " << tStepEnd << " time steps\n");

                        wavefields->resetWavefields();

                        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
                        Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
                        sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);

                        if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0) {
                            sources.getSeismogramHandler().filter(freqFilter);
                            receiversTrue.getSeismogramHandler().filter(freqFilter);
                        }

                        if (config.get<bool>("useSourceSignalInversion")) {
                            sourceEst.applyFilter(sources, shotNumber);
                            if (config.get<bool>("useSourceSignalTaper"))
                                sourceSignalTaper.apply(sources.getSeismogramHandler());
                        }

                        if (config.get<IndexType>("useSeismogramTaper") == 2) {                                                   
                            seismogramTaper2D.apply(receiversTrue.getSeismogramHandler()); 
                        }
                        seismogramTaper1D.apply(receiversTrue.getSeismogramHandler());

                        start_t_shot = common::Walltime::get();

                        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                            solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);
                        }

                        // check wavefield and seismogram for NaNs or infinite values
                        if ((commShot->any(!wavefields->isFinite(dist)) || commShot->any(!receivers.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write model and break
                            model->write("model_crash", config.get<IndexType>("FileFormat"));
                            COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield, output model as model_crash.FILE_EXTENSION!");
                        }
                        
//                        /* Normalize observed and synthetic data */
//                        if (config.get<bool>("NormalizeTraces")){
//                            receivers.getSeismogramHandler().normalize();
//                            receiversTrue.getSeismogramHandler().normalize();
//                        }

                        if (config.get<IndexType>("useSeismogramTaper") == 2) {                                                   
                            seismogramTaper2D.apply(receivers.getSeismogramHandler()); 
                        }
                        seismogramTaper1D.apply(receivers.getSeismogramHandler());
                            
                        receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);

                        /* Normalize observed and synthetic data */
                        if (config.get<bool>("NormalizeTraces")){
                            receivers.getSeismogramHandler().normalize();
                            receiversTrue.getSeismogramHandler().normalize();
                        }

                        /* Calculate misfit of one shot */
                        misfitPerIt.setValue(shotInd, dataMisfit->calc(receivers, receiversTrue));

                        end_t_shot = common::Walltime::get();
                        HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << numshots << ": Finished additional forward run\n");

                    } //end of loop over shots

                    commInterShot->sumArray(misfitPerIt.getLocalValues());
                    dataMisfit->addToStorage(misfitPerIt);

                    SLsearch.appendToLogFile(commAll, workflow.workflowStage + 1, workflow.iteration + 1, logFilename, dataMisfit->getMisfitSum(workflow.iteration + 1));
                    dataMisfit->clearStorage();
                    
                    if (useStreamConfig == 0) {
                        if (workflow.workflowStage != workflow.maxStage - 1) {
                            HOST_PRINT(commAll, "\nChange workflow stage\n");
                            workflow.changeStage(config, *dataMisfit, steplengthInit);
                        }
                    }

                } // end extra forward modelling

            } // end of loop over iterations
            
            if (useStreamConfig) {
                IndexType mincount = *std::min_element(std::begin(filterHistoryCount), std::end(filterHistoryCount));
                if (mincount == maxcount) {
                    HOST_PRINT(commAll, "\nMaximum number of iterations per pershot reached \n");
                    break;
                }
            }
            
        } // end of loop over pershot shots
        std::ofstream outFile("filterHistory.stage_" + std::to_string(workflow.workflowStage + 1) + ".txt");
        for (const auto &e : filterHistoryCount) outFile << e << "\n";
        
        if (useStreamConfig) {
            if (workflow.workflowStage != workflow.maxStage - 1) {
                HOST_PRINT(commAll, "\nChange workflow stage\n");
                workflow.changeStage(config, *dataMisfit, steplengthInit);
            }
        }

    } // end of loop over workflow stages

    
    globalEnd_t = common::Walltime::get();
    HOST_PRINT(commAll, "\nTotal runtime of WAVE-Inversion: " << globalEnd_t - globalStart_t << " sec.\nWAVE-Inversion finished!\n\n");
    return 0;
}
