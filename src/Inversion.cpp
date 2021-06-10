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

#include <AcquisitionEM/Receivers.hpp>
#include <AcquisitionEM/Sources.hpp>
#include <ForwardSolverEM/ForwardSolver.hpp>
#include <ForwardSolverEM/ForwardSolverFactory.hpp>
#include <ForwardSolverEM/SourceReceiverImpl/SourceReceiverImpl.hpp>
#include <ModelparameterEM/ModelparameterFactory.hpp>
#include <WavefieldsEM/WavefieldsFactory.hpp>

#include "Misfit/AbortCriterion.hpp"
#include "Misfit/Misfit.hpp"
#include "Misfit/MisfitFactory.hpp"
#include "Optimization/OptimizationFactory.hpp"
#include "Preconditioning/EnergyPreconditioning.hpp"
#include "SourceEstimation/SourceEstimation.hpp"
#include "StepLengthSearch/StepLengthSearch.hpp"
#include "Taper/Taper1D.hpp"
#include "Taper/Taper2D.hpp"

#include "Gradient/GradientCalculation.hpp"
#include "Gradient/GradientFactory.hpp"
#include "Workflow/Workflow.hpp"

#include "GradientEM/GradientCalculation.hpp"
#include "GradientEM/GradientFactory.hpp"
#include "WorkflowEM/Workflow.hpp"

#include <Common/HostPrint.hpp>

using namespace scai;
using namespace KITGPI;

bool verbose; // global variable definition

int main(int argc, char *argv[])
{
    double start_t, end_t, start_t_shot, end_t_shot; /* For timing */
    double globalStart_t, globalEnd_t; /* For timing */
    globalStart_t = common::Walltime::get();
    
    /* --------------------------------------- */
    /* Read configuration from file            */
    /* --------------------------------------- */
    Configuration::Configuration config;
    Configuration::Configuration configEM; 
    if (argc == 2) {
        config.readFromFile(argv[1], true);
        configEM.readFromFile(argv[1], true);
    } else if (argc == 3) {
        config.readFromFile(argv[1], true);
        configEM.readFromFile(argv[2], true);
    } else {
        std::cout << "\n\nNo configuration file given!\n\n" << std::endl;
        return (2);
    }
    verbose = config.get<bool>("verbose");
    
    // Initialization of inversion
    IndexType inversionType = config.get<IndexType>("inversionType");
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);   
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);  
    bool isSeismic = Common::checkEquationType<ValueType>(equationType);
    
    IndexType inversionTypeEM = configEM.get<IndexType>("inversionType");
    std::string dimensionEM = configEM.get<std::string>("dimension");
    std::string equationTypeEM = configEM.get<std::string>("equationType"); 
    std::transform(dimensionEM.begin(), dimensionEM.end(), dimensionEM.begin(), ::tolower);   
    std::transform(equationTypeEM.begin(), equationTypeEM.end(), equationTypeEM.begin(), ::tolower);  
    bool isSeismicEM = Common::checkEquationType<ValueType>(equationTypeEM);      
    
    // In case of that there is only one input configuration file
    if (argc == 2) {
        if (isSeismic) {
            inversionTypeEM = 0;
        } else {
            inversionType = 0;
        }
    }
    
    Configuration::Configuration configBig;
    Configuration::Configuration configBigEM;
    bool useStreamConfig = config.getAndCatch("useStreamConfig", false);
    bool useStreamConfigEM = configEM.getAndCatch("useStreamConfig", false);
    
    if (inversionType != 0 && useStreamConfig) {
        configBig.readFromFile(config.get<std::string>("streamConfigFilename"), true);
    }
    if (inversionTypeEM != 0 && useStreamConfigEM) {
        configBigEM.readFromFile(configEM.get<std::string>("streamConfigFilename"), true);    
    }
    
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
    IndexType tStepEndEM = static_cast<IndexType>((configEM.get<ValueType>("T") / configEM.get<ValueType>("DT")) + 0.5);  

    /* exchangeStrategy: 0 = self-restraint, 1 = self-restraint + mutual restraint (single parameter), 2 = self-restraint + mutual restraint (all parameters), 3 = self-restraint + mutual restraint (single parameter) + staged exchange, 4 = self-restraint + mutual restraint (all parameters) + sequential exchange, 5 = mutual restraint (single parameter) + staged exchange, 6 = mutual restraint (all parameters) + sequential exchange. */
    IndexType exchangeStrategy = config.get<IndexType>("exchangeStrategy");
    IndexType exchangeStrategyEM = configEM.get<IndexType>("exchangeStrategy");
    SCAI_ASSERT_ERROR(exchangeStrategy == exchangeStrategyEM, "exchangeStrategy != exchangeStrategyEM"); // check whether exchangeStrategy has been applied sucessfully.
        
    /*breakLoop: 0 = break iteration loop of which the abort criterion is true while do not influence another one. 0 is the same as individual inversion. 1 = break iteration loop if one of the two abort criterion is true. 2 = break iteration loop if both the two abort criterion is true. */
    IndexType breakLoopType = config.get<IndexType>("breakLoopType"); 
    IndexType breakLoopTypeEM = configEM.get<IndexType>("breakLoopType");
    SCAI_ASSERT_ERROR(breakLoopType == breakLoopTypeEM, "breakLoopType != breakLoopTypeEM"); // check whether breakLoopType has been applied sucessfully.
    if (inversionType == 1 || inversionTypeEM == 1 || exchangeStrategy > 2) {
        SCAI_ASSERT_ERROR(breakLoopType == 0, "breakLoopType != 0"); // check whether breakLoopType has been applied sucessfully.
    }

    std::string misfitType = config.get<std::string>("misfitType");
    std::string gradname(config.get<std::string>("gradientFilename"));
    std::string logFilename = config.get<std::string>("logFilename");
    ValueType steplengthInit = config.get<ValueType>("steplengthInit");
    std::string optimizationType = config.get<std::string>("optimizationType");
    
    std::string misfitTypeEM = configEM.get<std::string>("misfitType");
    std::string gradnameEM(configEM.get<std::string>("gradientFilename"));
    std::string logFilenameEM = configEM.get<std::string>("logFilename");
    ValueType steplengthInitEM = configEM.get<ValueType>("steplengthInit");
    std::string optimizationTypeEM = configEM.get<std::string>("optimizationType");
    
    IndexType maxiterations = config.get<IndexType>("maxIterations");
    if (inversionType == 0 && inversionTypeEM != 0) {
        maxiterations = configEM.get<IndexType>("maxIterations");
    } 
    
    /* inter node communicator */
    dmemo::CommunicatorPtr commAll = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(commAll->getNodeRank());
    
    if (inversionType != 0) {
        HOST_PRINT(commAll, "\n WAVE-Inversion " << dimension << " " << equationType << " - LAMA Version\n");
        HOST_PRINT(commAll, "","  - Running on " << commAll->getSize() << " mpi processes -\n\n");    
        if (commAll->getRank() == MASTERGPI) {
            config.print();
        }
    }
    if (inversionTypeEM != 0) {
        HOST_PRINT(commAll, "\n WAVE-Inversion " << dimensionEM << " " << equationTypeEM << " - LAMA Version\n");
        HOST_PRINT(commAll, "","  - Running on " << commAll->getSize() << " mpi processes -\n\n");    
        if (commAll->getRank() == MASTERGPI) {
            configEM.print();
        }
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
    
    Acquisition::Coordinates<ValueType> modelCoordinatesEM(configEM);
    Acquisition::Coordinates<ValueType> modelCoordinatesInversionEM(configEM, configEM.get<IndexType>("DHInversion"));
    
    Acquisition::Coordinates<ValueType> modelCoordinatesBig;
    Acquisition::Coordinates<ValueType> modelCoordinatesBigEM;

    if (inversionType != 0 && config.get<bool>("useVariableGrid")) {
        CheckParameter::checkVariableGrid(config, commAll, modelCoordinates);
        for (int layer=0;layer<modelCoordinates.getNumLayers();layer++){
            HOST_PRINT(commAll, "\n Number of gridpoints in layer: " << layer << " = " << modelCoordinates.getNGridpoints(layer));
        }
        auto numGridpointsRegular=config.get<IndexType>("NX")*config.get<IndexType>("NY")*config.get<IndexType>("NZ");
        HOST_PRINT(commAll, "\n Number of gripoints total: " << modelCoordinates.getNGridpoints());
        HOST_PRINT(commAll, "\n Percentage of gridpoints of the underlying regular grid given by NX*NY*NZ: "  << (float) modelCoordinates.getNGridpoints()/numGridpointsRegular * 100 << "% \n\n");
    }
    if (inversionTypeEM != 0 && configEM.get<bool>("useVariableGrid")) {
        CheckParameter::checkVariableGrid(configEM, commAll, modelCoordinatesEM);
        for (int layer=0;layer<modelCoordinatesEM.getNumLayers();layer++){
            HOST_PRINT(commAll, "\n Number of gridpoints in layer: " << layer << " = " << modelCoordinatesEM.getNGridpoints(layer)); 
        }
        auto numGridpointsRegular=configEM.get<IndexType>("NX")*configEM.get<IndexType>("NY")*configEM.get<IndexType>("NZ");
        HOST_PRINT(commAll, "\n Number of gridpoints total: " << modelCoordinatesEM.getNGridpoints());
        HOST_PRINT(commAll, "\n Percentage of gridpoints of the underlying regular grid given by NX*NY*NZ: "  << (float) modelCoordinatesEM.getNGridpoints()/numGridpointsRegular * 100 << "% \n\n");
    }

    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */    
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    IndexType shotDomain; // will contain the domain to which this processor belongs
    if (inversionType != 0) {
        shotDomain = Partitioning::getShotDomain(config, commAll);
    } else {
        shotDomain = Partitioning::getShotDomain(configEM, commAll); 
    }

    // Build subsets of processors for the shots
    dmemo::CommunicatorPtr commShot = commAll->split(shotDomain);
    dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());
    SCAI_DMEMO_TASK(commShot)

    /* --------------------------------------- */
    /* Distribution                            */
    /* --------------------------------------- */    
    // distribute the grid onto available processors
    dmemo::DistributionPtr dist = nullptr;
    dmemo::DistributionPtr distInversion = nullptr;
    dmemo::DistributionPtr distEM = nullptr;
    dmemo::DistributionPtr distInversionEM = nullptr;
    dmemo::DistributionPtr distBig = nullptr;
    dmemo::DistributionPtr distBigEM = nullptr;
    if (inversionType != 0) {
        if (config.get<IndexType>("partitioning") == 0 || config.get<IndexType>("partitioning") == 2) {
            //Block distribution = starting distribution for graph partitioner
            dist = std::make_shared<dmemo::BlockDistribution>(modelCoordinates.getNGridpoints(), commShot);
        } else if (config.get<IndexType>("partitioning") == 1) {
            SCAI_ASSERT(!config.get<bool>("useVariableGrid"), "Grid distribution is not available for the variable grid");
            dist = Partitioning::gridPartition<ValueType>(config, commShot);
            distInversion = Partitioning::gridPartitionInversion<ValueType>(config, commShot);
        } else {
            COMMON_THROWEXCEPTION("unknown partitioning method");
        }
        if (config.get<bool>("writeCoordinate") && shotDomain == 0) {
            // every shotdomain owns the same coordinates
            modelCoordinates.writeCoordinates(dist, ctx, config.get<std::string>("coordinateFilename"), config.get<IndexType>("FileFormat"));
        }
        if (useStreamConfig) {
            modelCoordinatesBig.init(configBig);
            distBig = Partitioning::gridPartition<ValueType>(configBig, commShot);
        }
    }
    if (inversionTypeEM != 0) {
        if (configEM.get<IndexType>("partitioning") == 0 || configEM.get<IndexType>("partitioning") == 2) {
            //Block distribution = starting distribution for graph partitioner
            distEM = std::make_shared<dmemo::BlockDistribution>(modelCoordinatesEM.getNGridpoints(), commShot);
        } else if (configEM.get<IndexType>("partitioning") == 1) {
            SCAI_ASSERT(!configEM.get<bool>("useVariableGrid"), "Grid distribution is not available for the variable grid");
            distEM = Partitioning::gridPartition<ValueType>(configEM, commShot);
            distInversionEM = Partitioning::gridPartitionInversion<ValueType>(configEM, commShot);
        } else {
            COMMON_THROWEXCEPTION("unknown partitioning method");
        }
        if (configEM.get<bool>("writeCoordinate") && shotDomain == 0) {
            // every shotdomain owns the same coordinates
            modelCoordinatesEM.writeCoordinates(distEM, ctx, configEM.get<std::string>("coordinateFilename"), configEM.get<IndexType>("FileFormat"));
        }
        if (useStreamConfigEM) {
            modelCoordinatesBigEM.init(configBigEM);
            distBigEM = Partitioning::gridPartition<ValueType>(configBigEM, commShot);
        }
    }
    
    /* --------------------------------------- */
    /* Factories                               */
    /* --------------------------------------- */
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivativesEM(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimensionEM));
    
    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));
    
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr modelPriori(Modelparameter::Factory<ValueType>::Create(equationType));
    
    typedef typename Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
    wavefieldPtr wavefields(Wavefields::Factory<ValueType>::Create(dimension, equationType));

    ForwardSolver::ForwardSolverEM<ValueType>::ForwardSolverPtr solverEM(ForwardSolver::FactoryEM<ValueType>::Create(dimensionEM, equationTypeEM));
    
    Modelparameter::ModelparameterEM<ValueType>::ModelparameterPtr modelEM(Modelparameter::FactoryEM<ValueType>::Create(equationTypeEM));
    Modelparameter::ModelparameterEM<ValueType>::ModelparameterPtr modelPrioriEM(Modelparameter::FactoryEM<ValueType>::Create(equationTypeEM));
    
    typedef typename Wavefields::WavefieldsEM<ValueType>::WavefieldPtr wavefieldPtrEM;
    wavefieldPtrEM wavefieldsEM(Wavefields::FactoryEM<ValueType>::Create(dimensionEM, equationTypeEM));
    
    /* --------------------------------------- */
    /* Memory estimation                       */
    /* --------------------------------------- */
    IndexType dtinversion = config.get<IndexType>("DTInversion");
    IndexType dtinversionEM = configEM.get<IndexType>("DTInversion");
    IndexType numShotDomains = config.get<IndexType>("NumShotDomains"); // total number of shot domains
    IndexType numShotDomainsEM = configEM.get<IndexType>("NumShotDomains");
    if (inversionType != 0) {
        ValueType memDerivatives = derivatives->estimateMemory(config, dist, modelCoordinates);
        ValueType memWavefileds = wavefields->estimateMemory(dist);
        ValueType memWavefiledsStorage;
        if (dimension.compare("3d") == 0) {
            memWavefiledsStorage = memWavefileds* (int) (tStepEnd / dtinversion) / pow(config.get<IndexType>("DHInversion"), 3);
        } else {
            memWavefiledsStorage = memWavefileds* (int) (tStepEnd / dtinversion) / pow(config.get<IndexType>("DHInversion"), 2);
        }
        ValueType memModel = model->estimateMemory(dist);
        ValueType memSolver = solver->estimateMemory(config, dist, modelCoordinates);
        ValueType memTotal = memDerivatives + memWavefileds + memModel + memSolver + memWavefiledsStorage;

        HOST_PRINT(commAll, "============== Memory Estimation " << equationType << ": ===============\n\n")
        HOST_PRINT(commAll, " -  Derivative Matrices \t" << memDerivatives << " MB\n");
        HOST_PRINT(commAll, " -  Wavefield vectors \t\t" << memWavefileds << " MB\n");
        HOST_PRINT(commAll, " -  Forward wavefield storage \t" << memWavefiledsStorage << " MB\n");
        HOST_PRINT(commAll, " -  Model Vectors \t\t" << memModel << " MB\n");
        HOST_PRINT(commAll, " -  Boundary Condition Vectors \t" << memSolver << " MB\n");
        HOST_PRINT(commAll, "\n Memory Usage (total / per partition): \n " << memTotal << " / " << memTotal / dist->getNumPartitions() << " MB ");
        IndexType numShotDomains = config.get<IndexType>("NumShotDomains"); // total number of shot domains
        if (numShotDomains > 1)
            HOST_PRINT(commAll, "\n Total Memory Usage (" << numShotDomains << " shot Domains ): \n " << memTotal * numShotDomains << " MB  ");

        HOST_PRINT(commAll, "\n\n ===================================================\n\n")    
    }    
    if (inversionTypeEM != 0) {
        ValueType memDerivativesEM = derivativesEM->estimateMemory(configEM, distEM, modelCoordinatesEM);
        ValueType memWavefiledsEM = wavefieldsEM->estimateMemory(distEM);
        ValueType memWavefiledsStorageEM;
        if (dimensionEM.compare("3d") == 0) {
            memWavefiledsStorageEM = memWavefiledsEM* (int) (tStepEndEM / dtinversionEM) / pow(configEM.get<IndexType>("DHInversion"), 3);
        } else {
            memWavefiledsStorageEM = memWavefiledsEM* (int) (tStepEndEM / dtinversionEM) / pow(configEM.get<IndexType>("DHInversion"), 2);
        }
        ValueType memModelEM = modelEM->estimateMemory(distEM);
        ValueType memSolverEM = solverEM->estimateMemory(configEM, distEM, modelCoordinatesEM);   
        ValueType memTotalEM = memDerivativesEM + memWavefiledsEM + memModelEM + memSolverEM + memWavefiledsStorageEM;

        HOST_PRINT(commAll, " =============== Memory Estimation " << equationTypeEM << ": =============\n\n");
        HOST_PRINT(commAll, "\n -  Derivative Matrices \t" << memDerivativesEM << " MB\n");
        HOST_PRINT(commAll, " -  Wavefield vectors \t\t" << memWavefiledsEM << " MB\n");
        HOST_PRINT(commAll, " -  Forward wavefield storage \t" << memWavefiledsStorageEM << " MB\n");
        HOST_PRINT(commAll, " -  Model Vectors \t\t" << memModelEM << " MB\n");
        HOST_PRINT(commAll, " -  Boundary Condition Vectors \t" << memSolverEM << " MB\n");
        HOST_PRINT(commAll, "\n Memory Usage (total / per partition): \n " << memTotalEM << " / " << memTotalEM / distEM->getNumPartitions() << " MB ");    
        
        if (numShotDomainsEM > 1) {      
            HOST_PRINT(commAll, "\n Total Memory Usage (" << numShotDomainsEM << " shot Domains ): \n " << memTotalEM * numShotDomainsEM << " MB  ");
        }        
        HOST_PRINT(commAll, "\n\n ===================================================\n\n")    
    }
    
    /* --------------------------------------- */
    /* Call partitioner                        */
    /* --------------------------------------- */
    if (inversionType != 0 && config.get<IndexType>("partitioning") == 2) {
        start_t = common::Walltime::get();
        dist = Partitioning::graphPartition(config, ctx, commShot, dist, *derivatives, modelCoordinates);
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "", "Finished graph partitioning in " << end_t - start_t << " sec.\n\n");
    }    
    
    if (inversionTypeEM != 0 && configEM.get<IndexType>("partitioning") == 2) {
        start_t = common::Walltime::get();
        distEM = Partitioning::graphPartition(configEM, ctx, commShot, distEM, *derivativesEM, modelCoordinatesEM);
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "", "Finished graph partitioning in " << end_t - start_t << " sec.\n\n");
    }    

    /* --------------------------------------- */
    /* Calculate derivative matrices           */
    /* --------------------------------------- */
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivativesInversionEM(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimensionEM));
    
    if (inversionType != 0) {
        start_t = common::Walltime::get();
        derivatives->init(dist, ctx, config, modelCoordinates, commShot);
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "", "Finished initializing matrices " << equationType << " in " << end_t - start_t << " sec.\n\n");
    }
    
    if (inversionTypeEM != 0) {
        start_t = common::Walltime::get();
        derivativesEM->init(distEM, ctx, configEM, modelCoordinatesEM, commShot);
        if (!useStreamConfigEM) {
            derivativesInversionEM->init(distEM, ctx, configEM, modelCoordinatesEM, commShot);
        } else {
            derivativesInversionEM->init(distBigEM, ctx, configEM, modelCoordinatesBigEM, commShot);
        }
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "", "Finished initializing matrices " << equationTypeEM << " in " << end_t - start_t << " sec.\n\n");
    }
            
    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */    
    if (inversionType != 0) {
        wavefields->init(ctx, dist);
    }    
    if (inversionTypeEM != 0) {
        wavefieldsEM->init(ctx, distEM);
    }

    /* --------------------------------------- */
    /* Wavefield record                        */
    /* --------------------------------------- */
    std::vector<wavefieldPtr> wavefieldrecord;
    std::vector<wavefieldPtrEM> wavefieldrecordEM;
    
    if (inversionType != 0) {
        for (IndexType i = 0; i < tStepEnd; i++) {
            if (i % dtinversion == 0) {
                // Only the local shared_ptr can be used to initialize a std::vector
                wavefieldPtr wavefieldsInversion = Wavefields::Factory<ValueType>::Create(dimension, equationType);
                wavefieldsInversion->init(ctx, distInversion);
                wavefieldrecord.push_back(wavefieldsInversion);
            }
        }
    }    
    if (inversionTypeEM != 0) {
        for (IndexType i = 0; i < tStepEndEM; i++) {
            if (i % dtinversionEM == 0) {
                wavefieldPtrEM wavefieldsInversionEM = Wavefields::FactoryEM<ValueType>::Create(dimensionEM, equationTypeEM);
                wavefieldsInversionEM->init(ctx, distInversionEM);
                wavefieldrecordEM.push_back(wavefieldsInversionEM);
            }
        }
    }
    
    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    Acquisition::Sources<ValueType> sources;
    Acquisition::Receivers<ValueType> receivers;
    Acquisition::SourcesEM<ValueType> sourcesEM;   
    Acquisition::ReceiversEM<ValueType> receiversEM;
    
    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings;
    std::vector<Acquisition::coordinate3D> cutCoordinates;
    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsEM;
    std::vector<Acquisition::coordinate3D> cutCoordinatesEM; 
    
    std::vector<scai::IndexType> uniqueShotNos;
    IndexType numshots = 1;
    IndexType numCuts = 1;
    std::vector<scai::IndexType> uniqueShotNosEM;
    IndexType numshotsEM = 1;
    IndexType numCutsEM = 1;
    
    std::shared_ptr<const dmemo::BlockDistribution> shotDist;
    std::shared_ptr<const dmemo::BlockDistribution> shotDistEM;
    IndexType maxcount = 1;   
    IndexType maxcountEM = 1;  
    
    if (inversionType != 0) {
        if (useStreamConfig) {
            std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;
            sources.getAcquisitionSettings(configBig, sourceSettingsBig);
            Acquisition::getCutCoord(cutCoordinates, sourceSettingsBig);
            Acquisition::getSettingsPerShot<ValueType>(sourceSettings, sourceSettingsBig, cutCoordinates);
        } else {
            sources.getAcquisitionSettings(config, sourceSettings);
        }
        HOST_PRINT(commAll, "\n================ checkSources " << equationType << " ===============\n");
        CheckParameter::checkSources(sourceSettings, modelCoordinates, commShot);
    
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
        numshots = uniqueShotNos.size();
        if (useStreamConfig) {
            numCuts = cutCoordinates.size();
            SCAI_ASSERT_ERROR(numshots == numCuts, "numshots != numCuts"); // check whether mdel pershot has been applied sucessfully.
            Acquisition::writeCutCoordToFile(config.get<std::string>("cutCoordinatesFilename"), cutCoordinates, uniqueShotNos);        
        }
        if (config.getAndCatch("useRandomSource", 0) != 0) {  
            shotDist = dmemo::blockDistribution(numShotDomains, commInterShot);
        } else {
            shotDist = dmemo::blockDistribution(numshots, commInterShot);
        }
        if (config.getAndCatch("useRandomSource", 0) != 0) { 
            maxcount = ceil((ValueType)maxiterations * numShotDomains / numshots * 1.2);
        }
        if (config.get<IndexType>("useReceiversPerShot") == 0) {
            receivers.init(config, modelCoordinates, ctx, dist);
        }
    }
    if (inversionTypeEM != 0) {
        if (useStreamConfigEM) {
            std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBigEM;
            sourcesEM.getAcquisitionSettings(configBigEM, sourceSettingsBigEM);
            Acquisition::getCutCoord(cutCoordinatesEM, sourceSettingsBigEM);
            Acquisition::getSettingsPerShot<ValueType>(sourceSettingsEM, sourceSettingsBigEM, cutCoordinatesEM);
        } else {
            sourcesEM.getAcquisitionSettings(configEM, sourceSettingsEM);
        }
        HOST_PRINT(commAll, "\n================ checkSources " << equationTypeEM << " ===============\n");
        CheckParameter::checkSources(sourceSettingsEM, modelCoordinatesEM, commShot);
    
        Acquisition::calcuniqueShotNo(uniqueShotNosEM, sourceSettingsEM);
        numshotsEM = uniqueShotNosEM.size();
        if (useStreamConfigEM) {
            numCutsEM = cutCoordinatesEM.size();
            SCAI_ASSERT_ERROR(numshotsEM == numCutsEM, "numshotsEM != numCutsEM"); // check whether mdel pershot has been applied sucessfully.
            Acquisition::writeCutCoordToFile(configEM.get<std::string>("cutCoordinatesFilename"), cutCoordinatesEM, uniqueShotNosEM);        
        }    
        if (configEM.getAndCatch("useRandomSource", 0) != 0) {  
            shotDistEM = dmemo::blockDistribution(numShotDomainsEM, commInterShot);
        } else {
            shotDistEM = dmemo::blockDistribution(numshotsEM, commInterShot);
        }     
        if (configEM.getAndCatch("useRandomSource", 0) != 0) { 
            maxcountEM = ceil((ValueType)maxiterations * numShotDomainsEM / numshotsEM * 1.2);
        }    
        if (configEM.get<IndexType>("useReceiversPerShot") == 0) {
            receiversEM.init(configEM, modelCoordinatesEM, ctx, distEM);
        }
    }
    
    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr modelPerShot(Modelparameter::Factory<ValueType>::Create(equationType));
    Modelparameter::ModelparameterEM<ValueType>::ModelparameterPtr modelPerShotEM(Modelparameter::FactoryEM<ValueType>::Create(equationTypeEM));
         
    if (inversionType != 0) {
        // load starting model
        if (!useStreamConfig) {    
            model->prepareForInversion(config, commShot);
            modelPriori->prepareForInversion(config, commShot);
            model->init(config, ctx, dist, modelCoordinates);
            modelPriori->init(config, ctx, dist, modelCoordinates);
        } else { 
            model->prepareForInversion(config, commShot);
            modelPriori->prepareForInversion(config, commShot);
            modelPerShot->prepareForInversion(config, commShot); // prepareForInversion is necessary for modelPerShot to calculate gradient.
            model->init(configBig, ctx, distBig, modelCoordinatesBig);   
            modelPriori->init(configBig, ctx, distBig, modelCoordinatesBig);
            modelPerShot->init(config, ctx, dist, modelCoordinates); 
        }
    }    
    if (inversionTypeEM != 0 || exchangeStrategy == 4 || exchangeStrategy == 6) {
        if (!useStreamConfigEM) {    
            modelEM->prepareForInversion(configEM, commShot);
            modelPrioriEM->prepareForInversion(configEM, commShot);
            modelEM->init(configEM, ctx, distEM, modelCoordinatesEM);
            modelPrioriEM->init(configEM, ctx, distEM, modelCoordinatesEM);
        } else {
            modelEM->prepareForInversion(configEM, commShot);
            modelPrioriEM->prepareForInversion(configEM, commShot);
            modelPerShotEM->prepareForInversion(configEM, commShot);// prepareForInversion is necessary for modelPerShot to calculate gradient.
            modelEM->init(configBigEM, ctx, distBigEM, modelCoordinatesBigEM);  
            modelPrioriEM->init(configBigEM, ctx, distBigEM, modelCoordinatesBigEM);
            modelPerShotEM->init(configEM, ctx, distEM, modelCoordinatesEM);  
        }
    }

    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */
    if (inversionType != 0) {
        if (!useStreamConfig) {
            solver->initForwardSolver(config, *derivatives, *wavefields, *model, modelCoordinates, ctx, config.get<ValueType>("DT"));
        } else {
            solver->initForwardSolver(config, *derivatives, *wavefields, *modelPerShot, modelCoordinates, ctx, config.get<ValueType>("DT"));
        }
    }
    if (inversionTypeEM != 0) {
        if (!useStreamConfigEM) {
            solverEM->initForwardSolver(configEM, *derivativesEM, *wavefieldsEM, *modelEM, modelCoordinatesEM, ctx, configEM.get<ValueType>("DT"));
        } else {
            solverEM->initForwardSolver(configEM, *derivativesEM, *wavefieldsEM, *modelPerShotEM, modelCoordinatesEM, ctx, configEM.get<ValueType>("DT"));
        }
    }

    /* --------------------------------------- */
    /* True data                               */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> receiversTrue;
    Acquisition::ReceiversEM<ValueType> receiversTrueEM;
    Acquisition::Receivers<ValueType> receiversStart;
    Acquisition::ReceiversEM<ValueType> receiversStartEM;
    Taper::Taper2D<ValueType> seismogramTaper2D;
    Taper::Taper1D<ValueType> seismogramTaper1D;
    Taper::Taper2D<ValueType> seismogramTaper2DEM;
    Taper::Taper1D<ValueType> seismogramTaper1DEM;
    
    if (inversionType != 0) {
        if (config.get<IndexType>("useReceiversPerShot") == 0) {
            receiversTrue.init(config, modelCoordinates, ctx, dist);
            receiversStart.init(config, modelCoordinates, ctx, dist);
        }
        seismogramTaper1D.init(std::make_shared<dmemo::NoDistribution>(tStepEnd), ctx, 1);
    }
    if (inversionTypeEM != 0) {
        if (configEM.get<IndexType>("useReceiversPerShot") == 0) {
            receiversTrueEM.init(configEM, modelCoordinatesEM, ctx, distEM);
            receiversStartEM.init(configEM, modelCoordinatesEM, ctx, distEM);
        }
        seismogramTaper1DEM.init(std::make_shared<dmemo::NoDistribution>(tStepEndEM), ctx, 1);
    }
    
    /* --------------------------------------- */
    /* Misfit                                  */
    /* --------------------------------------- */
    Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(Misfit::Factory<ValueType>::Create(misfitType));
    lama::DenseVector<ValueType> misfitPerIt(numshots, 0, ctx);
    ValueType weightingStabilizingFunctionalGradient = 1.0;
    ValueType weightingCrossGradient = 1.0;
    
    Misfit::Misfit<ValueType>::MisfitPtr dataMisfitEM(Misfit::Factory<ValueType>::Create(misfitTypeEM));
    lama::DenseVector<ValueType> misfitPerItEM(numshotsEM, 0, ctx);  
    ValueType weightingStabilizingFunctionalGradientEM = 1.0;   
    ValueType weightingCrossGradientEM = 1.0;

    /* --------------------------------------- */
    /* Adjoint sources                         */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> adjointSources;
    Acquisition::ReceiversEM<ValueType> adjointSourcesEM;
    if (inversionType != 0 && config.get<IndexType>("useReceiversPerShot") == 0)
        adjointSources.init(config, modelCoordinates, ctx, dist);

    if (inversionTypeEM != 0 && configEM.get<IndexType>("useReceiversPerShot") == 0)
        adjointSourcesEM.init(configEM, modelCoordinatesEM, ctx, distEM);
    
    /* --------------------------------------- */
    /* Workflow                                */
    /* --------------------------------------- */
    Workflow::Workflow<ValueType> workflow(config);
    Workflow::WorkflowEM<ValueType> workflowEM(configEM);

    /* --------------------------------------- */
    /* Abort criterion                         */
    /* --------------------------------------- */
    AbortCriterion<ValueType> abortCriterion;
    bool breakLoop = false;
    AbortCriterion<ValueType> abortCriterionEM;
    bool breakLoopEM = false;

    /* --------------------------------------- */
    /* Step length search                      */
    /* --------------------------------------- */
    StepLengthSearch<ValueType> SLsearch;
    StepLengthSearch<ValueType> SLsearchEM;
    if (inversionType != 0)
        SLsearch.initLogFile(commAll, logFilename, misfitType);
    if (inversionTypeEM != 0)
        SLsearchEM.initLogFile(commAll, logFilenameEM, misfitTypeEM);

    /* --------------------------------------- */
    /* Source estimation                       */
    /* --------------------------------------- */
    SourceEstimation<ValueType> sourceEst;
    Taper::Taper1D<ValueType> sourceSignalTaper;
    SourceEstimation<ValueType> sourceEstEM;
    Taper::Taper1D<ValueType> sourceSignalTaperEM;
    
    if (inversionType != 0) {
        // calculate source dist
        scai::lama::DenseVector<IndexType> sourcecoords = getsourcecoordinates(sourceSettings, modelCoordinates);
        scai::dmemo::DistributionPtr dist_sources = Acquisition::calcDistribution(sourcecoords, dist);
        if (config.get<bool>("useSourceSignalInversion"))
            sourceEst.init(config, ctx, dist_sources, sourceSignalTaper);
    }
    if (inversionTypeEM != 0) {
        // calculate source dist
        scai::lama::DenseVector<IndexType> sourcecoordsEM = getsourcecoordinates(sourceSettingsEM, modelCoordinatesEM);
        scai::dmemo::DistributionPtr dist_sourcesEM = Acquisition::calcDistribution(sourcecoordsEM, distEM);
        if (configEM.get<bool>("useSourceSignalInversion"))
            sourceEstEM.init(configEM, ctx, dist_sourcesEM, sourceSignalTaperEM);
    }
    
    /* --------------------------------------- */
    /* Frequency filter                        */
    /* --------------------------------------- */
    Filter::Filter<ValueType> freqFilter;
    std::string transFcnFmly = "butterworth";
    if (inversionType != 0 && (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0))
        freqFilter.init(config.get<ValueType>("DT"), tStepEnd);

    Filter::Filter<ValueType> freqFilterEM;
    if (inversionTypeEM != 0 && (workflowEM.getLowerCornerFreq() != 0.0 || workflowEM.getUpperCornerFreq() != 0.0))
        freqFilterEM.init(configEM.get<ValueType>("DT"), tStepEndEM);
    
    /* --------------------------------------- */
    /* Gradients                               */
    /* --------------------------------------- */
    typedef typename Gradient::Gradient<ValueType>::GradientPtr gradientPtr;
    gradientPtr gradient(Gradient::Factory<ValueType>::Create(equationType));
    gradientPtr gradientPerShot(Gradient::Factory<ValueType>::Create(equationType));
    gradientPtr stabilizingFunctionalGradient(Gradient::Factory<ValueType>::Create(equationType));
    gradientPtr crossGradientDerivative(Gradient::Factory<ValueType>::Create(equationType));
    
    typedef typename Gradient::GradientEM<ValueType>::GradientPtr gradientPtrEM;
    gradientPtrEM gradientEM(Gradient::FactoryEM<ValueType>::Create(equationTypeEM));
    gradientPtrEM gradientPerShotEM(Gradient::FactoryEM<ValueType>::Create(equationTypeEM));
    gradientPtrEM stabilizingFunctionalGradientEM(Gradient::FactoryEM<ValueType>::Create(equationTypeEM));
    gradientPtrEM crossGradientDerivativeEM(Gradient::FactoryEM<ValueType>::Create(equationTypeEM)); 
            
    if (inversionType != 0) {
        if (!useStreamConfig) {
            gradient->init(ctx, dist);
            gradientPerShot->init(ctx, dist);       
            stabilizingFunctionalGradient->init(ctx, dist);
            crossGradientDerivative->init(ctx, dist);  
        } else {
            gradient->init(ctx, distBig);
            stabilizingFunctionalGradient->init(ctx, distBig);
            crossGradientDerivative->init(ctx, distBig);
        }
        gradient->setNormalizeGradient(config.get<bool>("normalizeGradient"));
        gradientPerShot->setNormalizeGradient(config.get<bool>("normalizeGradient"));  
        stabilizingFunctionalGradient->setNormalizeGradient(config.get<bool>("normalizeGradient"));
        crossGradientDerivative->setNormalizeGradient(config.get<bool>("normalizeGradient"));     
    }
    if (inversionTypeEM != 0) {
        if (!useStreamConfigEM) {
            gradientEM->init(ctx, distEM);
            gradientPerShotEM->init(ctx, distEM);      
            stabilizingFunctionalGradientEM->init(ctx, distEM); 
            crossGradientDerivativeEM->init(ctx, distEM);  
        } else {
            gradientEM->init(ctx, distBigEM);
            stabilizingFunctionalGradientEM->init(ctx, distBigEM);
            crossGradientDerivativeEM->init(ctx, distBigEM);
        }
        gradientEM->setNormalizeGradient(configEM.get<bool>("normalizeGradient"));
        gradientPerShotEM->setNormalizeGradient(configEM.get<bool>("normalizeGradient"));  
        stabilizingFunctionalGradientEM->setNormalizeGradient(configEM.get<bool>("normalizeGradient")); 
        crossGradientDerivativeEM->setNormalizeGradient(configEM.get<bool>("normalizeGradient"));     
    }
    
    /* --------------------------------------- */
    /* Gradient calculation                    */
    /* --------------------------------------- */
    GradientCalculation<ValueType> gradientCalculation;
    GradientCalculationEM<ValueType> gradientCalculationEM;

    /* --------------------------------------- */
    /* Gradient taper                          */
    /* --------------------------------------- */
    Preconditioning::SourceReceiverTaper<ValueType> ReceiverTaper;
    Taper::Taper1D<ValueType> gradientTaper1D;
    Taper::Taper2D<ValueType> wavefieldTaper2D;
    Preconditioning::SourceReceiverTaper<ValueType> ReceiverTaperEM;
    Taper::Taper1D<ValueType> gradientTaper1DEM; 
    Taper::Taper2D<ValueType> wavefieldTaper2DEM;
    Taper::Taper2D<ValueType> modelTaper2DJoint;
    
    if (inversionType != 0) {
        if (config.get<IndexType>("useReceiversPerShot") == 0)
            ReceiverTaper.init(dist, ctx, receivers, config, modelCoordinates, config.get<IndexType>("receiverTaperRadius"));
        if (config.get<bool>("useGradientTaper")) {
            if (!useStreamConfig) {
                gradientTaper1D.init(dist, ctx, 1);
                gradientTaper1D.read(config.get<std::string>("gradientTaperName"), config.get<IndexType>("FileFormat"));
            } else {
                gradientTaper1D.init(distBig, ctx, 1);
                gradientTaper1D.read(configBig.get<std::string>("gradientTaperName"), config.get<IndexType>("FileFormat"));
            }  
        }
        wavefieldTaper2D.initWavefieldTransform(config, distInversion, dist, ctx, isSeismic); 
        wavefieldTaper2D.calcInversionAverageMatrix(modelCoordinates, modelCoordinatesInversion);  
    }
    if (inversionTypeEM != 0) {
        if (configEM.get<IndexType>("useReceiversPerShot") == 0)
            ReceiverTaperEM.init(distEM, ctx, receiversEM, configEM, modelCoordinatesEM, configEM.get<IndexType>("receiverTaperRadius"));
        if (configEM.get<bool>("useGradientTaper")) {
            if (!useStreamConfigEM) {   
                gradientTaper1DEM.init(distEM, ctx, 1);
                gradientTaper1DEM.read(configEM.get<std::string>("gradientTaperName"), configEM.get<IndexType>("FileFormat"));
            } else {        
                gradientTaper1DEM.init(distBigEM, ctx, 1);
                gradientTaper1DEM.read(configBigEM.get<std::string>("gradientTaperName"), configEM.get<IndexType>("FileFormat"));
            }
        }
        wavefieldTaper2DEM.initWavefieldTransform(configEM, distInversionEM, distEM, ctx, isSeismicEM);
        wavefieldTaper2DEM.calcInversionAverageMatrix(modelCoordinatesEM, modelCoordinatesInversionEM);
    } 
    if (inversionType != 0 && inversionTypeEM != 0) {
        if (!useStreamConfig) {
            modelTaper2DJoint.initModelTransform(dist, distEM, ctx);  
            modelTaper2DJoint.calcSeismictoEMMatrix(modelCoordinates, config, modelCoordinatesEM, configEM);
            modelTaper2DJoint.calcEMtoSeismicMatrix(modelCoordinates, config, modelCoordinatesEM, configEM); 
        } else {
            modelTaper2DJoint.initModelTransform(distBig, distBigEM, ctx);  
            modelTaper2DJoint.calcSeismictoEMMatrix(modelCoordinatesBig, configBig, modelCoordinatesBigEM, configBigEM);
            modelTaper2DJoint.calcEMtoSeismicMatrix(modelCoordinatesBig, configBig, modelCoordinatesBigEM, configBigEM); 
        }
    }
    
    /* --------------------------------------- */
    /* Gradient preconditioning                */
    /* --------------------------------------- */
    Preconditioning::EnergyPreconditioning<ValueType> energyPrecond;
    if (inversionType != 0) {
        energyPrecond.init(dist, config, isSeismic);
    }
    Preconditioning::EnergyPreconditioning<ValueType> energyPrecondEM;
    if (inversionTypeEM != 0) {
        energyPrecondEM.init(distEM, configEM, isSeismicEM);
    }
    
    /* --------------------------------------- */
    /* Gradient optimization                   */
    /* --------------------------------------- */
    Optimization::Optimization<ValueType>::OptimizationPtr gradientOptimization(Optimization::Factory<ValueType>::Create(optimizationType));
    Optimization::Optimization<ValueType>::OptimizationPtr gradientOptimizationEM(Optimization::Factory<ValueType>::Create(optimizationTypeEM));
    
    if (inversionType != 0) 
        gradientOptimization->init(dist);
    if (inversionTypeEM != 0) 
        gradientOptimizationEM->init(distEM);  
    
    /* --------------------------------------- */
    /*       Modelparameter preparation        */
    /* --------------------------------------- */    
    if (inversionType != 0 && (inversionType == 3 || model->getParameterisation() == 2 || model->getParameterisation() == 1)) {
        // in case of using porosity and saturation in inversion
        model->calcRockMatrixParameter(config);     
    }
    if (inversionTypeEM != 0 && (inversionTypeEM == 3 || modelEM->getParameterisation() == 2 || modelEM->getParameterisation() == 1)) {
        // in case of using porosity and saturation in inversion
        modelEM->calcRockMatrixParameter(configEM);     
    }
    
    /* --------------------------------------- */
    /*       Loop over workflow stages         */
    /* --------------------------------------- */
    for (workflow.workflowStage = 0; workflow.workflowStage < workflow.maxStage; workflow.workflowStage++) {

        if (inversionType != 0) {        
            breakLoop = false;
            workflow.printParameters(commAll);

            gradientCalculation.allocate(config, dist, ctx, workflow);
            seismogramTaper1D.calcTimeDampingTaper(workflow.getTimeDampingFactor(), config.get<ValueType>("DT"));  

            if (workflow.getLowerCornerFreq() != 0.0 && workflow.getUpperCornerFreq() != 0.0)
                freqFilter.calc(transFcnFmly, "bp", workflow.getFilterOrder(), workflow.getLowerCornerFreq(), workflow.getUpperCornerFreq());
            else if (workflow.getLowerCornerFreq() != 0.0 && workflow.getUpperCornerFreq() == 0.0)
                freqFilter.calc(transFcnFmly, "lp", workflow.getFilterOrder(), workflow.getLowerCornerFreq());
            else if (workflow.getLowerCornerFreq() == 0.0 && workflow.getUpperCornerFreq() != 0.0)
                freqFilter.calc(transFcnFmly, "hp", workflow.getFilterOrder(), workflow.getUpperCornerFreq());            
        }
        
        if (inversionTypeEM != 0) {
            breakLoopEM = false;
            workflowEM.workflowStage = workflow.workflowStage;
            workflowEM.printParameters(commAll);

            gradientCalculationEM.allocate(configEM, distEM, ctx, workflowEM);
            seismogramTaper1DEM.calcTimeDampingTaper(workflowEM.getTimeDampingFactor(), configEM.get<ValueType>("DT"));

            if (workflowEM.getLowerCornerFreq() != 0.0 && workflowEM.getUpperCornerFreq() != 0.0)
                freqFilterEM.calc(transFcnFmly, "bp", workflowEM.getFilterOrder(), workflowEM.getLowerCornerFreq(), workflowEM.getUpperCornerFreq());
            else if (workflowEM.getLowerCornerFreq() != 0.0 && workflowEM.getUpperCornerFreq() == 0.0)
                freqFilterEM.calc(transFcnFmly, "lp", workflowEM.getFilterOrder(), workflowEM.getLowerCornerFreq());
            else if (workflowEM.getLowerCornerFreq() == 0.0 && workflowEM.getUpperCornerFreq() != 0.0)
                freqFilterEM.calc(transFcnFmly, "hp", workflowEM.getFilterOrder(), workflowEM.getUpperCornerFreq());                                
        }
        
        /* --------------------------------------- */
        /*        Loop over iterations             */
        /* --------------------------------------- */ 
        std::vector<scai::IndexType> uniqueShotNosRand(numShotDomains, 0); 
        std::vector<scai::IndexType> uniqueShotNosRandEM(numShotDomainsEM, 0);
        std::vector<scai::IndexType> shotHistory(numshots, 0);
        std::vector<scai::IndexType> shotHistoryEM(numshotsEM, 0);
        std::vector<scai::IndexType> misfitTypeHistory(misfitType.length() - 2, 0);
        std::vector<scai::IndexType> misfitTypeHistoryEM(misfitTypeEM.length() - 2, 0);
        IndexType shotNumber;
        IndexType shotIndTrue = 0;            
        for (workflow.iteration = 0; workflow.iteration < maxiterations; workflow.iteration++) {

            /* --------------------------------------- */
            /*        Start the first inversion        */
            /* --------------------------------------- */            
            if (inversionType != 0 && (breakLoop == false || breakLoopType == 2)) {
                // Begin of one seismic model update 
                HOST_PRINT(commAll, "\n=================================================");
                HOST_PRINT(commAll, "\n=========== Workflow stage " << workflow.workflowStage + 1 << " of " << workflow.maxStage << " ===============");
                HOST_PRINT(commAll, "\n============     Iteration " << workflow.iteration + 1 << "       ==============");
                HOST_PRINT(commAll, "\n=================== " << equationType << " =======================\n\n");
                start_t = common::Walltime::get();
                
                /* Update model for fd simulation (averaging, inverse Density ...) */
                if (!useStreamConfig) { 
                    model->prepareForModelling(modelCoordinates, ctx, dist, commShot);  
                    solver->prepareForModelling(*model, config.get<ValueType>("DT"));
                }                 
                
                if (workflow.iteration == 0 && commInterShot->getRank() == 0) {
                    /* only shot Domain 0 writes output */
                    model->write((config.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration)), config.get<IndexType>("FileFormat"));
                }
                
                /* --------------------------------------- */
                /*        Loop over Seismic shots          */
                /* --------------------------------------- */
                gradient->resetGradient(); // reset gradient because gradient is a sum of all gradientsPerShot gradients+=gradientPerShot
                crossGradientDerivative->resetGradient();
                misfitPerIt = 0;

                if (config.getAndCatch("useRandomSource", 0) != 0) { 
                    start_t = common::Walltime::get();
                    Acquisition::getRandShotNos<ValueType>(uniqueShotNosRand, shotHistory, uniqueShotNos, maxcount, config.getAndCatch("useRandomSource", 0));
                    Acquisition::writeRandShotNosToFile(commAll, logFilename, uniqueShotNosRand, workflow.workflowStage + 1, workflow.iteration + 1, config.getAndCatch("useRandomSource", 0));
                    end_t = common::Walltime::get();
                    HOST_PRINT(commAll, "Finished initializing a random shot sequence: " << workflow.iteration + 1 << " of " << maxiterations << " (maxcount = " << maxcount << ") in " << end_t - start_t << " sec.\n");
                }
                dataMisfit->init(config, misfitTypeHistory, numshots); // in case of that random misfit function is used
                IndexType localShotInd = 0;     
                for (IndexType shotInd = shotDist->lb(); shotInd < shotDist->ub(); shotInd++) {
                    if (config.getAndCatch("useRandomSource", 0) == 0) {  
                        shotNumber = uniqueShotNos[shotInd];
                        shotIndTrue = shotInd;
                    } else {
                        shotNumber = uniqueShotNosRand[shotInd];
                        Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNos, shotNumber);
                    }
                    localShotInd++;
                    HOST_PRINT(commShot, "Shot number " << shotNumber << ", local shot " << localShotInd << " of " << shotDist->getLocalSize() << ": started\n");
                        
                    if (useStreamConfig) {
                        HOST_PRINT(commShot, "Switch to model shot: " << shotIndTrue + 1 << " of " << numshots << "\n");
                        model->getModelPerShot(*modelPerShot, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotIndTrue)); 
                        modelPerShot->prepareForModelling(modelCoordinates, ctx, dist, commShot); 
                        solver->prepareForModelling(*modelPerShot, config.get<ValueType>("DT"));
                    }

                    if (config.get<IndexType>("useReceiversPerShot") == 1) {
                        receivers.init(config, modelCoordinates, ctx, dist, shotNumber);
                        receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber);
                        receiversStart.init(config, modelCoordinates, ctx, dist, shotNumber);
                        adjointSources.init(config, modelCoordinates, ctx, dist, shotNumber);

                        ReceiverTaper.init(dist, ctx, receivers, config, modelCoordinates, config.get<IndexType>("receiverTaperRadius"));
                    } else if (config.get<IndexType>("useReceiversPerShot") == 2) {
                        receivers.init(config, modelCoordinates, ctx, dist, shotNumber, numshots);
                        receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber, numshots);
                        receiversStart.init(config, modelCoordinates, ctx, dist, shotNumber, numshots);
                        adjointSources.init(config, modelCoordinates, ctx, dist, shotNumber, numshots);

                        ReceiverTaper.init(dist, ctx, receivers, config, modelCoordinates, config.get<IndexType>("receiverTaperRadius"));
                    }

                    /* Read field data (or pseudo-observed data, respectively) */
                    receiversTrue.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".shot_" + std::to_string(shotNumber), 1);
                    receiversStart.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber), 1);
                                        
                    if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0){
                        receiversTrue.getSeismogramHandler().filter(freqFilter);
                        receiversStart.getSeismogramHandler().filter(freqFilter);
                    }

                    /* Reset approximated Hessian per shot */
                    energyPrecond.resetApproxHessian();

                    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
                    Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
                    sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);

                    if (!useStreamConfig) {
                        CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(config, sourceSettingsShot, *model, modelCoordinates, shotNumber);
                    } else {
                        CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(config, sourceSettingsShot, *modelPerShot, modelCoordinates, shotNumber);
                    }
                    
                    if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0)
                        sources.getSeismogramHandler().filter(freqFilter);
                    
                    /* Source time function inversion */
                    if (config.get<bool>("useSourceSignalInversion")){
                        if (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1) {
                            HOST_PRINT(commShot, "Shot number " << shotNumber << ", local shot " << localShotInd << " of " << shotDist->getLocalSize() << " : Source Time Function Inversion\n");

                            wavefields->resetWavefields();

                            if (!useStreamConfig) {
                                for (scai::IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                                    solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);
                                }
                            } else {
                                for (scai::IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                                    solver->run(receivers, sources, *modelPerShot, *wavefields, *derivatives, tStep);
                                }
                            }

                            solver->resetCPML();
                            
                            /* Normalize observed and synthetic data */
                            if (config.get<bool>("normalizeTraces")){
                                receivers.getSeismogramHandler().normalize();
                                receiversTrue.getSeismogramHandler().normalize();
                            }

                            if (config.get<bool>("maxOffsetSrcEst"))
                                sourceEst.calcOffsetMutes(sources, receivers, config.get<ValueType>("maxOffsetSrcEst"),modelCoordinates);
                            
                            sourceEst.estimateSourceSignal(receivers, receiversTrue, shotIndTrue, shotNumber);

                            sourceEst.applyFilter(sources, shotIndTrue);
                            if (config.get<bool>("useSourceSignalTaper"))
                                sourceSignalTaper.apply(sources.getSeismogramHandler());
                            if (config.get<bool>("writeInvertedSource"))
                                sources.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("sourceSeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
                        } else {
                            sourceEst.applyFilter(sources, shotIndTrue);
                            if (config.get<bool>("useSourceSignalTaper"))
                                sourceSignalTaper.apply(sources.getSeismogramHandler());
                        }
                    }
                    
                    if (config.get<IndexType>("useSeismogramTaper") > 1) {
                        seismogramTaper2D.init(receiversTrue.getSeismogramHandler());
                        if (config.get<IndexType>("useSeismogramTaper") == 4) {
                            seismogramTaper2D.read(config.get<std::string>("seismogramTaperName") + ".misfitCalc.shot_" + std::to_string(shotNumber) + ".mtx");
                        } else {
                            seismogramTaper2D.read(config.get<std::string>("seismogramTaperName") + ".shot_" + std::to_string(shotNumber) + ".mtx");
                        }
                        seismogramTaper2D.apply(receiversTrue.getSeismogramHandler()); 
                        seismogramTaper2D.apply(receiversStart.getSeismogramHandler()); 
                    }
                    seismogramTaper1D.apply(receiversTrue.getSeismogramHandler());
                    seismogramTaper1D.apply(receiversStart.getSeismogramHandler());
                    
                    if (workflow.iteration == 0){
                        receiversTrue.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
                        receiversStart.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
                    }
                    
                    /* --------------------------------------- */
                    /*        Forward modelling 1              */
                    /* --------------------------------------- */
                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshots << ": Start time stepping with " << tStepEnd << " time steps\n");

                    wavefields->resetWavefields();

                    start_t_shot = common::Walltime::get();
                    
                    if (!useStreamConfig) {
                        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {

                            solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);
                            if (tStep % dtinversion == 0) {
                                // save wavefields in std::vector
                                *wavefieldrecord[floor(tStep / dtinversion + 0.5)] = wavefieldTaper2D.applyWavefieldAverage(wavefields);
                            }
                        }
                    } else {
                        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {

                            solver->run(receivers, sources, *modelPerShot, *wavefields, *derivatives, tStep);
                            if (tStep % dtinversion == 0) {
                                // save wavefields in std::vector
                                *wavefieldrecord[floor(tStep / dtinversion + 0.5)] = wavefieldTaper2D.applyWavefieldAverage(wavefields);
                            }
                        }
                    }

                    // check wavefield and seismogram for NaNs or infinite values
                    if ((commShot->any(!wavefields->isFinite(dist)) || commShot->any(!receivers.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write model and break
                        model->write("model_crash", config.get<IndexType>("FileFormat"));
                    COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield, output model as model_crash.FILE_EXTENSION!");
                    }

                    if (config.get<IndexType>("useSeismogramTaper") == 2 || config.get<IndexType>("useSeismogramTaper") == 3) {                                                   
                        seismogramTaper2D.apply(receivers.getSeismogramHandler()); 
                    }
                    seismogramTaper1D.apply(receivers.getSeismogramHandler());
                                    
                    receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration) + ".shot_" + std::to_string(shotNumber), modelCoordinates);

                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshots << ": Calculate misfit and adjoint sources\n");

                    /* Normalize observed and synthetic data */
                    if (config.get<bool>("normalizeTraces")){
                        receivers.getSeismogramHandler().normalize();
                        receiversTrue.getSeismogramHandler().normalize();
                    }
                    
                    /* Calculate misfit of one shot */
                    misfitPerIt.setValue(shotIndTrue, dataMisfit->calc(receivers, receiversTrue, shotIndTrue));

                    /* Calculate adjoint sources */
                    dataMisfit->calcAdjointSources(adjointSources, receivers, receiversTrue, shotIndTrue);

                    /* Calculate gradient */
                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshots << ": Start Backward\n");
                    if (!useStreamConfig) {
                        gradientCalculation.run(commAll, *solver, *derivatives, receivers, sources, adjointSources, *model, *gradientPerShot, wavefieldrecord, config, modelCoordinates, shotNumber, workflow, wavefieldTaper2D, energyPrecond);
                    } else {
                        gradientCalculation.run(commAll, *solver, *derivatives, receivers, sources, adjointSources, *modelPerShot, *gradientPerShot, wavefieldrecord, config, modelCoordinates, shotNumber, workflow, wavefieldTaper2D, energyPrecond);
                    }

                    /* Apply energy preconditioning per shot */
                    energyPrecond.apply(*gradientPerShot, shotNumber, config.get<IndexType>("FileFormat"));
                    
                    if (config.get<IndexType>("useReceiversPerShot") != 0)
                        ReceiverTaper.apply(*gradientPerShot);

                    gradientPerShot->normalize();

                    if (!useStreamConfig) {
                        *gradient += *gradientPerShot;
                    } else {
                        gradient->sumGradientPerShot(*gradientPerShot, modelCoordinates, modelCoordinatesBig, cutCoordinates, shotIndTrue, config.get<IndexType>("BoundaryWidth"));
                    }
                    
                    solver->resetCPML();

                    end_t_shot = common::Walltime::get();
                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshots << ": Finished in " << end_t_shot - start_t_shot << " sec.\n");

                } //end of loop over shots

                gradient->sumShotDomain(commInterShot);                
                commInterShot->sumArray(misfitPerIt.getLocalValues());

                HOST_PRINT(commAll, "\n======== Finished loop over shots " << equationType << " =========");
                HOST_PRINT(commAll, "\n=================================================\n");

                dataMisfit->addToStorage(misfitPerIt);
                misfitPerIt = 0;
                if (misfitType.length() > 2) {
                    dataMisfit->appendMisfitTypeShotsToFile(commAll, logFilename, workflow.workflowStage + 1, workflow.iteration + 1);
                    dataMisfit->appendMisfitsToFile(commAll, logFilename, workflow.workflowStage + 1, workflow.iteration + 1);
                }

                SLsearch.appendToLogFile(commAll, workflow.workflowStage + 1, workflow.iteration, logFilename, dataMisfit->getMisfitSum(workflow.iteration));

                /* Check abort criteria */
                HOST_PRINT(commAll, "\nMisfit after stage " << workflow.workflowStage + 1 << ", iteration " << workflow.iteration << ": " << dataMisfit->getMisfitSum(workflow.iteration) << "\n");
            
                /* Apply receiver Taper (if ReceiverTaperRadius=0 gradient will be multplied by 1) */
                if (config.get<IndexType>("useReceiversPerShot") == 0)
                    ReceiverTaper.apply(*gradient);
                
                if (config.get<bool>("useGradientTaper"))
                    gradientTaper1D.apply(*gradient);

                gradient->applyMedianFilter(config);  
                
                if (inversionType == 2) {                        
                    
                    // joint inversion with cross-gradient constraint
                    HOST_PRINT(commAll, "\n===========================================");
                    HOST_PRINT(commAll, "\n========= calcCrossGradient " << equationType << " =========");
                    HOST_PRINT(commAll, "\n===========================================\n");
                    gradient->normalize();  
                    
                    crossGradientDerivativeEM->calcModelDerivative(*dataMisfitEM, *modelEM, *derivativesInversionEM, configEM, workflowEM);
                    
                    crossGradientDerivative->calcCrossGradient(*dataMisfitEM, *model, *derivativesInversionEM, configEM, modelTaper2DJoint, workflow);  
                    if (config.get<bool>("useGradientTaper")) {
                        gradientTaper1D.apply(*crossGradientDerivative);
                    }  
                    crossGradientDerivative->applyMedianFilter(config); 
                    crossGradientDerivative->normalize();  
                    if (config.get<IndexType>("writeGradient") == 2 && commInterShot->getRank() == 0){
                        crossGradientDerivative->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".crossGradient", config.get<IndexType>("FileFormat"), workflow);
                    }
                                            
                    crossGradientDerivative->calcCrossGradientDerivative(*dataMisfitEM, *model, *derivativesInversionEM, configEM, modelTaper2DJoint, workflow);
                    if (config.get<bool>("useGradientTaper")) {
                        gradientTaper1D.apply(*crossGradientDerivative);
                    }  
                    crossGradientDerivative->applyMedianFilter(config); 
                    crossGradientDerivative->normalize();              
                    if (config.get<IndexType>("writeGradient") == 2 && commInterShot->getRank() == 0){
                        crossGradientDerivative->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".crossGradientDerivative", config.get<IndexType>("FileFormat"), workflow);
                    }      
                    if ((workflow.iteration > 0) && (dataMisfit->getMisfitSum(workflow.iteration - 1) - dataMisfit->getMisfitSum(workflow.iteration) - 0.01 * dataMisfit->getMisfitSum(workflow.iteration - 1 ) < 0)) {
                        weightingCrossGradient *= 0.8;
                    }
                    HOST_PRINT(commAll, "weightingCrossGradient " << equationType << " = " << weightingCrossGradient << "\n");  
                    HOST_PRINT(commAll, "\n===========================================\n");
                    
                    *crossGradientDerivative *= weightingCrossGradient;            
                    *gradient += *crossGradientDerivative;
                }
                
                if (workflow.iteration != 0 && config.get<IndexType>("stablizingFunctionalType") != 0) {
                    
                    // inversion with regularization constraint
                    HOST_PRINT(commAll, "\n===========================================");
                    HOST_PRINT(commAll, "\n=== calcStabilizingFunctionalGradient " << equationType << " ===");
                    HOST_PRINT(commAll, "\n===========================================\n");
                    gradient->normalize();  
                    
                    stabilizingFunctionalGradient->calcStabilizingFunctionalGradient(*model, *modelPriori, config, *dataMisfit, workflow);
                    if (config.get<bool>("useGradientTaper")) {
                        gradientTaper1D.apply(*stabilizingFunctionalGradient);
                    }  
                    stabilizingFunctionalGradient->applyMedianFilter(config);  
                    stabilizingFunctionalGradient->normalize();  
                    if (config.get<IndexType>("writeGradient") == 2 && commInterShot->getRank() == 0){
                        stabilizingFunctionalGradient->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".stabilizingFunctionalGradient", config.get<IndexType>("FileFormat"), workflow);
                    }  
                    if ((workflow.iteration > 0) && (dataMisfit->getMisfitSum(workflow.iteration - 1) - dataMisfit->getMisfitSum(workflow.iteration) - 0.01 * dataMisfit->getMisfitSum(workflow.iteration - 1 ) < 0)) {
                        weightingStabilizingFunctionalGradient *= 0.6;
                    }
                    HOST_PRINT(commAll, "weightingStabilizingFunctionalGradient = " << weightingStabilizingFunctionalGradient << "\n");  
                    HOST_PRINT(commAll, "\n===========================================\n");
                    
                    *stabilizingFunctionalGradient *= weightingStabilizingFunctionalGradient;
                    *gradient += *stabilizingFunctionalGradient;   
                }
                
                // scale function in gradientOptimization must be the final operation for gradient.
                gradientOptimization->apply(*gradient, workflow, *model, config);
                
                /* Output of gradient */
                /* only shot Domain 0 writes output */
                if (config.get<IndexType>("writeGradient") == 1 && commInterShot->getRank() == 0) {
                    gradient->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1), config.get<IndexType>("FileFormat"), workflow);
                }
                
                /* --------------------------------------- */
                /* Check abort criteria for two inversions */
                /* --------------------------------------- */              
                HOST_PRINT(commAll, "\n========== Check abort criteria " << equationType << " ==========\n"); 
                
                if (config.getAndCatch("useRandomSource", 0) == 0 && misfitType.length() == 2) { 
                    breakLoop = abortCriterion.check(commAll, *dataMisfit, config, steplengthInit, workflow, breakLoopEM, breakLoopType);
                }
                // We set a new break condition so that two inversions can change stage simutanously in joint inversion
                if (breakLoop == true) {                 
                    if (inversionType == 3 && exchangeStrategy != 0 && exchangeStrategy != 4 && exchangeStrategy != 6) {   
                        HOST_PRINT(commAll, "\n=================================================");
                        HOST_PRINT(commAll, "\n========= Joint petrophysical inversion =========");
                        HOST_PRINT(commAll, "\n=============  From Seismic to EM  ==============");
                        HOST_PRINT(commAll, "\n=================================================\n");
                        modelTaper2DJoint.exchangePetrophysics(*model, *modelEM, configEM);
                    }       
                    if (breakLoopType == 0 && breakLoopEM == true) {
                        break;
                    }           
                } 

                if (breakLoop == false || breakLoopType == 2) {
                    HOST_PRINT(commAll, "\n================================================");
                    HOST_PRINT(commAll, "\n========== Start step length search " << equationType << " ======\n");
                    SLsearch.run(commAll, *solver, *derivatives, receivers, sourceSettings, receiversTrue, *model, dist, config, modelCoordinates, *gradient, steplengthInit, dataMisfit->getMisfitIt(workflow.iteration), workflow, freqFilter, sourceEst, sourceSignalTaper, uniqueShotNosRand, dataMisfit->getMisfitTypeShots());

                    HOST_PRINT(commAll, "================= Update Model " << equationType << " ============\n\n");
                    /* Apply model update */
                    *gradient *= SLsearch.getSteplength();
                    *model -= *gradient;
                    
                    if (config.get<bool>("useModelThresholds"))
                        model->applyThresholds(config);

                    if (model->getParameterisation() == 2 || model->getParameterisation() == 1) {
                        HOST_PRINT(commAll, "\n======= calcWaveModulusFromPetrophysics " << equationType << " =======\n");  
                        model->calcWaveModulusFromPetrophysics();  
                        if (config.get<bool>("useModelThresholds"))
                            model->applyThresholds(config); 
                    } else if (inversionType == 3) {
                        HOST_PRINT(commAll, "\n======= calcPetrophysicsFromWaveModulus " << equationType << " =======\n");  
                        model->calcPetrophysicsFromWaveModulus();
                        if (config.get<bool>("useModelThresholds"))
                            model->applyThresholds(config); 
                        if (exchangeStrategy != 5 && exchangeStrategy != 6) {
                            // case 0,1,2,3,4: self-constraint of the petrophysical relationship    
                            HOST_PRINT(commAll, "\n======= calcWaveModulusFromPetrophysics " << equationType << " =======\n");  
                            model->calcWaveModulusFromPetrophysics(); 
                            if (config.get<bool>("useModelThresholds"))
                                model->applyThresholds(config); 
                        }                    
                    }
                    if (commInterShot->getRank() == 0) {
                        /* only shot Domain 0 writes output */
                        model->write((config.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1)), config.get<IndexType>("FileFormat"));
                    }
                    
                    steplengthInit *= 0.98;

                    end_t = common::Walltime::get();
                    HOST_PRINT(commAll, "\nFinished iteration " << workflow.iteration + 1 << " in " << end_t - start_t << " sec.\n\n");
                }
                
            }  // End of once Seismic gradient, dataMisfit calculation and model update

            if (inversionType == 3 && (exchangeStrategy == 1 || exchangeStrategy == 2 || (workflow.iteration == maxiterations - 1 && (exchangeStrategy == 3 || exchangeStrategy == 5)))) {
                HOST_PRINT(commAll, "\n=================================================");
                HOST_PRINT(commAll, "\n========= Joint petrophysical inversion =========");
                HOST_PRINT(commAll, "\n=============  From Seismic to EM  ==============");
                HOST_PRINT(commAll, "\n=================================================\n");                
                modelTaper2DJoint.exchangePetrophysics(*model, *modelEM, configEM);
            }
            
            /* --------------------------------------- */
            /*        Start the second inversion       */
            /* --------------------------------------- */            
            if (inversionTypeEM != 0 && (breakLoopEM == false || breakLoopType == 2)) {
                // Begin of one EM model update 
                workflowEM.iteration = workflow.iteration;
                HOST_PRINT(commAll, "\n=================================================");
                HOST_PRINT(commAll, "\n=========== Workflow stage " << workflowEM.workflowStage + 1 << " of " << workflowEM.maxStage << " ===============");
                HOST_PRINT(commAll, "\n============     Iteration " << workflowEM.iteration + 1 << "       ==============");
                HOST_PRINT(commAll, "\n=================== " << equationTypeEM << " =======================\n\n");
                start_t = common::Walltime::get();
                
                /* Update modelEM for fd simulation (averaging, getVelocityEM ...) */
                if (!useStreamConfigEM) { 
                    modelEM->prepareForModelling(modelCoordinatesEM, ctx, distEM, commShot);  
                    solverEM->prepareForModelling(*modelEM, configEM.get<ValueType>("DT"));
                }                 
                
                if (workflowEM.iteration == 0 && commInterShot->getRank() == 0) {
                    /* only shot Domain 0 writes output */
                    modelEM->write((configEM.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".It_" + std::to_string(workflowEM.iteration)), configEM.get<IndexType>("FileFormat"));
                }
                
                /* --------------------------------------- */
                /*        Loop over EM shots               */
                /* --------------------------------------- */
                gradientEM->resetGradient(); // reset gradient because gradient is a sum of all gradientsPerShot gradients+=gradientPerShot
                crossGradientDerivativeEM->resetGradient();
                misfitPerItEM = 0;

                if (configEM.getAndCatch("useRandomSource", 0) != 0) { 
                    start_t = common::Walltime::get();
                    Acquisition::getRandShotNos<ValueType>(uniqueShotNosRandEM, shotHistoryEM, uniqueShotNosEM, maxcountEM, configEM.getAndCatch("useRandomSource", 0));
                    Acquisition::writeRandShotNosToFile(commAll, logFilenameEM, uniqueShotNosRandEM, workflowEM.workflowStage + 1, workflowEM.iteration + 1, configEM.getAndCatch("useRandomSource", 0));
                    end_t = common::Walltime::get();
                    HOST_PRINT(commAll, "Finished initializing a random shot sequence: " << workflowEM.iteration + 1 << " of " << maxiterations << " (maxcount = " << maxcountEM << ") in " << end_t - start_t << " sec.\n");
                }
                dataMisfitEM->init(configEM, misfitTypeHistoryEM, numshotsEM); // in case of that random misfit function is used
                IndexType localShotInd = 0;     
                for (IndexType shotInd = shotDistEM->lb(); shotInd < shotDistEM->ub(); shotInd++) {
                    if (configEM.getAndCatch("useRandomSource", 0) == 0) {  
                        shotNumber = uniqueShotNosEM[shotInd];
                        shotIndTrue = shotInd;
                    } else {
                        shotNumber = uniqueShotNosRandEM[shotInd];
                        Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNosEM, shotNumber);
                    }
                    localShotInd++;
                    HOST_PRINT(commShot, "Shot number " << shotNumber << ", local shot " << localShotInd << " of " << shotDistEM->getLocalSize() << ": started\n");
                        
                    if (useStreamConfigEM) {
                        HOST_PRINT(commShot, "Switch to model shot: " << shotIndTrue + 1 << " of " << numshotsEM << "\n");
                        modelEM->getModelPerShot(*modelPerShotEM, modelCoordinatesEM, modelCoordinatesBigEM, cutCoordinatesEM.at(shotIndTrue)); 
                        modelPerShotEM->prepareForModelling(modelCoordinatesEM, ctx, distEM, commShot); 
                        solverEM->prepareForModelling(*modelPerShotEM, configEM.get<ValueType>("DT"));
                    }

                    if (configEM.get<IndexType>("useReceiversPerShot") == 1) {
                        receiversEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber);
                        receiversTrueEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber);
                        receiversStartEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber);
                        adjointSourcesEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber);

                        ReceiverTaperEM.init(distEM, ctx, receiversEM, configEM, modelCoordinatesEM, configEM.get<IndexType>("receiverTaperRadius"));
                    } else if (configEM.get<IndexType>("useReceiversPerShot") == 2) {
                        receiversEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber, numshotsEM);
                        receiversTrueEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber, numshotsEM);
                        receiversStartEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber, numshotsEM);
                        adjointSourcesEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber, numshotsEM);

                        ReceiverTaperEM.init(distEM, ctx, receiversEM, configEM, modelCoordinatesEM, configEM.get<IndexType>("receiverTaperRadius"));
                    }

                    /* Read field data (or pseudo-observed data, respectively) */
                    receiversTrueEM.getSeismogramHandler().read(configEM.get<IndexType>("SeismogramFormat"), configEM.get<std::string>("fieldSeisName") + ".shot_" + std::to_string(shotNumber), 1);
                    receiversStartEM.getSeismogramHandler().read(configEM.get<IndexType>("SeismogramFormat"), configEM.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber), 1);
                                        
                    if (workflowEM.getLowerCornerFreq() != 0.0 || workflowEM.getUpperCornerFreq() != 0.0){
                        receiversTrueEM.getSeismogramHandler().filter(freqFilterEM);
                        receiversStartEM.getSeismogramHandler().filter(freqFilterEM);
                    }

                    /* Reset approximated Hessian per shot */
                    energyPrecondEM.resetApproxHessian();

                    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
                    Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettingsEM, shotNumber);
                    sourcesEM.init(sourceSettingsShot, configEM, modelCoordinatesEM, ctx, distEM);

                    if (!useStreamConfigEM) {
                        CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(configEM, sourceSettingsShot, *modelEM, modelCoordinatesEM, shotNumber);
                    } else {
                        CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(configEM, sourceSettingsShot, *modelPerShotEM, modelCoordinatesEM, shotNumber);
                    }
                    
                    if (workflowEM.getLowerCornerFreq() != 0.0 || workflowEM.getUpperCornerFreq() != 0.0)
                        sourcesEM.getSeismogramHandler().filter(freqFilterEM);
                    
                    /* Source time function inversion */
                    if (configEM.get<bool>("useSourceSignalInversion")){
                        if (workflowEM.iteration == 0 || shotHistoryEM[shotIndTrue] == 1) {
                            HOST_PRINT(commShot, "Shot number " << shotNumber << ", local shot " << localShotInd << " of " << shotDistEM->getLocalSize() << " : Source Time Function Inversion\n");

                            wavefieldsEM->resetWavefields();

                            if (!useStreamConfigEM) {
                                for (scai::IndexType tStep = 0; tStep < tStepEndEM; tStep++) {
                                    solverEM->run(receiversEM, sourcesEM, *modelEM, *wavefieldsEM, *derivativesEM, tStep);
                                }
                            } else {
                                for (scai::IndexType tStep = 0; tStep < tStepEndEM; tStep++) {
                                    solverEM->run(receiversEM, sourcesEM, *modelPerShotEM, *wavefieldsEM, *derivativesEM, tStep);
                                }
                            }

                            solverEM->resetCPML();
                            
                            /* Normalize observed and synthetic data */
                            if (configEM.get<bool>("normalizeTraces")){
                                receiversEM.getSeismogramHandler().normalize();
                                receiversTrueEM.getSeismogramHandler().normalize();
                            }

                            if (configEM.get<bool>("maxOffsetSrcEst"))
                                sourceEstEM.calcOffsetMutes(sourcesEM, receiversEM, configEM.get<ValueType>("maxOffsetSrcEst"),modelCoordinatesEM);
                            
                            sourceEstEM.estimateSourceSignal(receiversEM, receiversTrueEM, shotIndTrue, shotNumber);

                            sourceEstEM.applyFilter(sourcesEM, shotIndTrue);
                            if (configEM.get<bool>("useSourceSignalTaper"))
                                sourceSignalTaperEM.apply(sourcesEM.getSeismogramHandler());
                            if (configEM.get<bool>("writeInvertedSource"))
                                sourcesEM.getSeismogramHandler().write(configEM.get<IndexType>("SeismogramFormat"), configEM.get<std::string>("sourceSeismogramFilename") + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinatesEM);
                        } else {
                            sourceEstEM.applyFilter(sourcesEM, shotIndTrue);
                            if (configEM.get<bool>("useSourceSignalTaper"))
                                sourceSignalTaperEM.apply(sourcesEM.getSeismogramHandler());
                        }
                    }
                    
                    if (configEM.get<IndexType>("useSeismogramTaper") > 1) {
                        seismogramTaper2DEM.init(receiversTrueEM.getSeismogramHandler());
                        if (configEM.get<IndexType>("useSeismogramTaper") == 4) {
                            seismogramTaper2DEM.read(configEM.get<std::string>("seismogramTaperName") + "misfitCalc.shot_" + std::to_string(shotNumber) + ".mtx");
                        } else {
                            seismogramTaper2DEM.read(configEM.get<std::string>("seismogramTaperName") + ".shot_" + std::to_string(shotNumber) + ".mtx");
                        }
                        seismogramTaper2DEM.apply(receiversTrueEM.getSeismogramHandler()); 
                        seismogramTaper2DEM.apply(receiversStartEM.getSeismogramHandler()); 
                    }
                    seismogramTaper1DEM.apply(receiversTrueEM.getSeismogramHandler());
                    seismogramTaper1DEM.apply(receiversStartEM.getSeismogramHandler());
                    
                    if (workflowEM.iteration == 0){
                        receiversTrueEM.getSeismogramHandler().write(configEM.get<IndexType>("SeismogramFormat"), configEM.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinatesEM);
                        receiversStartEM.getSeismogramHandler().write(configEM.get<IndexType>("SeismogramFormat"), configEM.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinatesEM);
                    }
                    
                    /* --------------------------------------- */
                    /*        Forward modelling 1              */
                    /* --------------------------------------- */
                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshotsEM << ": Start time stepping with " << tStepEndEM << " time steps\n");

                    wavefieldsEM->resetWavefields();

                    start_t_shot = common::Walltime::get();
                    
                    if (!useStreamConfigEM) {
                        for (IndexType tStep = 0; tStep < tStepEndEM; tStep++) {

                            solverEM->run(receiversEM, sourcesEM, *modelEM, *wavefieldsEM, *derivativesEM, tStep);
                            if (tStep % dtinversion == 0) {
                                // save wavefieldsEM in std::vector
                                *wavefieldrecordEM[floor(tStep / dtinversion + 0.5)] = wavefieldTaper2DEM.applyWavefieldAverage(wavefieldsEM);
                            }
                        }
                    } else {
                        for (IndexType tStep = 0; tStep < tStepEndEM; tStep++) {

                            solverEM->run(receiversEM, sourcesEM, *modelPerShotEM, *wavefieldsEM, *derivativesEM, tStep);
                            if (tStep % dtinversion == 0) {
                                // save wavefieldsEM in std::vector
                                *wavefieldrecordEM[floor(tStep / dtinversion + 0.5)] = wavefieldTaper2DEM.applyWavefieldAverage(wavefieldsEM);
                            }
                        }
                    }

                    // check wavefield and seismogram for NaNs or infinite values
                    if ((commShot->any(!wavefieldsEM->isFinite(distEM)) || commShot->any(!receiversEM.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write modelEM and break
                        modelEM->write("model_crash", configEM.get<IndexType>("FileFormat"));
                    COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield, output modelEM as model_crash.FILE_EXTENSION!");
                    }

                    if (configEM.get<IndexType>("useSeismogramTaper") == 2 || configEM.get<IndexType>("useSeismogramTaper") == 3) {                                                   
                        seismogramTaper2DEM.apply(receiversEM.getSeismogramHandler()); 
                    }
                    seismogramTaper1DEM.apply(receiversEM.getSeismogramHandler());
                                    
                    receiversEM.getSeismogramHandler().write(configEM.get<IndexType>("SeismogramFormat"), configEM.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".It_" + std::to_string(workflowEM.iteration) + ".shot_" + std::to_string(shotNumber), modelCoordinatesEM);

                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshotsEM << ": Calculate misfit and adjoint sources\n");

                    /* Normalize observed and synthetic data */
                    if (configEM.get<bool>("normalizeTraces")){
                        receiversEM.getSeismogramHandler().normalize();
                        receiversTrueEM.getSeismogramHandler().normalize();
                    }
                    
                    /* Calculate misfit of one shot */
                    misfitPerItEM.setValue(shotIndTrue, dataMisfitEM->calc(receiversEM, receiversTrueEM, shotIndTrue));

                    /* Calculate adjoint sources */
                    dataMisfitEM->calcAdjointSources(adjointSourcesEM, receiversEM, receiversTrueEM, shotIndTrue);

                    /* Calculate gradient */
                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshotsEM << ": Start Backward\n");
                    if (!useStreamConfigEM) {
                        gradientCalculationEM.run(commAll, *solverEM, *derivativesEM, receiversEM, sourcesEM, adjointSourcesEM, *modelEM, *gradientPerShotEM, wavefieldrecordEM, configEM, modelCoordinatesEM, shotNumber, workflowEM, wavefieldTaper2DEM, energyPrecondEM);
                    } else {
                        gradientCalculationEM.run(commAll, *solverEM, *derivativesEM, receiversEM, sourcesEM, adjointSourcesEM, *modelPerShotEM, *gradientPerShotEM, wavefieldrecordEM, configEM, modelCoordinatesEM, shotNumber, workflowEM, wavefieldTaper2DEM, energyPrecondEM);
                    }

                    /* Apply energy preconditioning per shot */
                    energyPrecondEM.apply(*gradientPerShotEM, shotNumber, configEM.get<IndexType>("FileFormat"));
                    
                    if (configEM.get<IndexType>("useReceiversPerShot") != 0)
                        ReceiverTaperEM.apply(*gradientPerShotEM);

                    gradientPerShotEM->normalize();

                    if (!useStreamConfigEM) {
                        *gradientEM += *gradientPerShotEM;
                    } else {
                        gradientEM->sumGradientPerShot(*gradientPerShotEM, modelCoordinatesEM, modelCoordinatesBigEM, cutCoordinatesEM, shotIndTrue, configEM.get<IndexType>("BoundaryWidth"));
                    }
                    
                    solverEM->resetCPML();

                    end_t_shot = common::Walltime::get();
                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshotsEM << ": Finished in " << end_t_shot - start_t_shot << " sec.\n");

                } //end of loop over shots

                dataMisfitEM->sumShotDomain(commInterShot);   
                gradientEM->sumShotDomain(commInterShot);                
                commInterShot->sumArray(misfitPerItEM.getLocalValues());

                HOST_PRINT(commAll, "\n======== Finished loop over shots " << equationTypeEM << " =========");
                HOST_PRINT(commAll, "\n=================================================\n");

                dataMisfitEM->addToStorage(misfitPerItEM);
                misfitPerItEM = 0;
                if (misfitTypeEM.length() > 2) {
                    dataMisfitEM->appendMisfitTypeShotsToFile(commAll, logFilenameEM, workflowEM.workflowStage + 1, workflowEM.iteration + 1);
                    dataMisfitEM->appendMisfitsToFile(commAll, logFilenameEM, workflowEM.workflowStage + 1, workflowEM.iteration + 1);
                }

                SLsearchEM.appendToLogFile(commAll, workflowEM.workflowStage + 1, workflowEM.iteration, logFilenameEM, dataMisfitEM->getMisfitSum(workflowEM.iteration));

                /* Check abort criteria */
                HOST_PRINT(commAll, "\nMisfit after stage " << workflowEM.workflowStage + 1 << ", iteration " << workflowEM.iteration << ": " << dataMisfitEM->getMisfitSum(workflowEM.iteration) << "\n");
            
                /* Apply receiver Taper (if ReceiverTaperRadius=0 gradient will be multplied by 1) */
                if (configEM.get<IndexType>("useReceiversPerShot") == 0)
                    ReceiverTaperEM.apply(*gradientEM);

                if (configEM.get<bool>("useGradientTaper"))
                    gradientTaper1DEM.apply(*gradientEM);
                
                gradientEM->applyMedianFilter(configEM);  

                if (inversionTypeEM == 2) {
                    
                    // joint inversion with cross-gradient constraint
                    HOST_PRINT(commAll, "\n===========================================");
                    HOST_PRINT(commAll, "\n========== calcCrossGradient " << equationTypeEM << " ========");
                    HOST_PRINT(commAll, "\n===========================================\n");
                    gradientEM->normalize();  
                    
                    crossGradientDerivative->calcModelDerivative(*dataMisfit, *model, *derivativesInversionEM, configEM, modelTaper2DJoint, workflow);
                    
                    crossGradientDerivativeEM->calcCrossGradient(*dataMisfit, *modelEM, *derivativesInversionEM, configEM, modelTaper2DJoint, workflowEM);
                    if (configEM.get<bool>("useGradientTaper")) {
                        gradientTaper1DEM.apply(*crossGradientDerivativeEM);
                    }  
                    crossGradientDerivativeEM->applyMedianFilter(configEM);  
                    crossGradientDerivativeEM->normalize();  
                    if (configEM.get<IndexType>("writeGradient") == 2 && commInterShot->getRank() == 0){
                        crossGradientDerivativeEM->write(gradnameEM + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".It_" + std::to_string(workflowEM.iteration + 1) + ".crossGradient", configEM.get<IndexType>("FileFormat"), workflowEM);
                    }  
                    
                    crossGradientDerivativeEM->calcCrossGradientDerivative(*dataMisfit, *modelEM, *derivativesInversionEM, configEM, modelTaper2DJoint, workflowEM); 
                    if (configEM.get<bool>("useGradientTaper")) {
                        gradientTaper1DEM.apply(*crossGradientDerivativeEM);
                    }   
                    crossGradientDerivativeEM->applyMedianFilter(configEM);      
                    crossGradientDerivativeEM->normalize();   
                    if (configEM.get<IndexType>("writeGradient") == 2 && commInterShot->getRank() == 0){
                        crossGradientDerivativeEM->write(gradnameEM + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".It_" + std::to_string(workflowEM.iteration + 1) + ".crossGradientDerivative", configEM.get<IndexType>("FileFormat"), workflowEM);
                    }             
                    if ((workflowEM.iteration > 0) && (dataMisfitEM->getMisfitSum(workflowEM.iteration - 1) - dataMisfitEM->getMisfitSum(workflowEM.iteration) - 0.01 * dataMisfitEM->getMisfitSum(workflowEM.iteration - 1 ) < 0)) {
                        weightingCrossGradientEM *= 0.8;
                    }
                    HOST_PRINT(commAll, "weightingCrossGradient " << equationTypeEM << " = " << weightingCrossGradientEM << "\n");  
                    HOST_PRINT(commAll, "\n===========================================\n");
                    
                    *crossGradientDerivativeEM *= weightingCrossGradientEM;
                    *gradientEM += *crossGradientDerivativeEM;           
                }
                
                if (workflowEM.iteration != 0 && configEM.get<IndexType>("stablizingFunctionalType") != 0) {
                    
                    // inversion with regularization constraint
                    HOST_PRINT(commAll, "\n===========================================");
                    HOST_PRINT(commAll, "\n=== calcStabilizingFunctionalGradient " << equationTypeEM << " ===");
                    HOST_PRINT(commAll, "\n===========================================\n");
                    gradientEM->normalize();  
                    
                    stabilizingFunctionalGradientEM->calcStabilizingFunctionalGradient(*modelEM, *modelPrioriEM, configEM, *dataMisfitEM, workflowEM);
                    if (config.get<bool>("useGradientTaper")) {
                        gradientTaper1DEM.apply(*stabilizingFunctionalGradientEM);
                    }  
                    stabilizingFunctionalGradientEM->applyMedianFilter(configEM);  
                    stabilizingFunctionalGradientEM->normalize();  
                    if (configEM.get<IndexType>("writeGradient") == 2 && commInterShot->getRank() == 0){
                        stabilizingFunctionalGradientEM->write(gradnameEM + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".It_" + std::to_string(workflowEM.iteration + 1) + ".stabilizingFunctionalGradientEM", configEM.get<IndexType>("FileFormat"), workflowEM);
                    }  
                    if ((workflowEM.iteration > 0) && (dataMisfitEM->getMisfitSum(workflowEM.iteration - 1) - dataMisfitEM->getMisfitSum(workflowEM.iteration) - 0.01 * dataMisfitEM->getMisfitSum(workflowEM.iteration - 1 ) < 0)) {
                        weightingStabilizingFunctionalGradientEM *= 0.6;
                    }
                    HOST_PRINT(commAll, "weightingStabilizingFunctionalGradient = " << weightingStabilizingFunctionalGradientEM << "\n");  
                    HOST_PRINT(commAll, "\n===========================================\n");
                    
                    *stabilizingFunctionalGradientEM *= weightingStabilizingFunctionalGradientEM;
                    *gradientEM += *stabilizingFunctionalGradientEM;   
                }
                
                gradientOptimizationEM->apply(*gradientEM, workflowEM, *modelEM, configEM);

                /* Output of gradient */
                /* only shot Domain 0 writes output */
                if (configEM.get<IndexType>("writeGradient") == 1 && commInterShot->getRank() == 0) {
                    gradientEM->write(gradnameEM + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".It_" + std::to_string(workflowEM.iteration + 1), configEM.get<IndexType>("FileFormat"), workflowEM);
                }
                
                /* --------------------------------------- */
                /* Check abort criteria for two inversions */
                /* --------------------------------------- */   
                HOST_PRINT(commAll, "\n========== Check abort criteria " << equationTypeEM << " ==========\n");  
                       
                if (configEM.getAndCatch("useRandomSource", 0) == 0 && misfitTypeEM.length() == 2) {     
                    breakLoopEM = abortCriterionEM.check(commAll, *dataMisfitEM, configEM, steplengthInitEM, workflowEM, breakLoop, breakLoopType);
                }
                // We set a new break condition so that two inversions can change stage simutanously in joint inversion
                if (breakLoopEM == true) {   
                    if (inversionTypeEM == 3 && exchangeStrategy != 0 && exchangeStrategy != 4 && exchangeStrategy != 6) {
                        HOST_PRINT(commAll, "\n=================================================");
                        HOST_PRINT(commAll, "\n========= Joint petrophysical inversion =========");
                        HOST_PRINT(commAll, "\n=============  From EM to Seismic  ==============");
                        HOST_PRINT(commAll, "\n=================================================\n");                    
                        modelTaper2DJoint.exchangePetrophysics(*modelEM, *model, config); 
                    }   
                    if (breakLoopType == 1) {
                        if (breakLoop == false) {
                            if(workflow.workflowStage != workflow.maxStage-1) {
                                HOST_PRINT(commAll, "\nChange workflow stage\n");
                                workflow.changeStage(config, *dataMisfit, steplengthInit);                
                            }
                        }
                        break;  
                    } else if (breakLoopType == 2) {
                        if (breakLoop == true) {
                            if(workflow.workflowStage != workflow.maxStage-1) {
                                HOST_PRINT(commAll, "\nChange workflow stage\n");
                                workflow.changeStage(config, *dataMisfit, steplengthInit);                
                            }
                            break; 
                        } else {
                            breakLoopEM = false;
                        }
                    } else if (breakLoopType == 0) {
                        if (breakLoop == true) {
                            break;
                        } else {
                            if (inversionType == 3 && exchangeStrategy != 0 && exchangeStrategy != 4 && exchangeStrategy != 6) {
                                HOST_PRINT(commAll, "\n=================================================");
                                HOST_PRINT(commAll, "\n========= Joint petrophysical inversion =========");
                                HOST_PRINT(commAll, "\n=============  From Seismic to EM  ==============");
                                HOST_PRINT(commAll, "\n=================================================\n");                        
                                modelTaper2DJoint.exchangePetrophysics(*model, *modelEM, configEM);
                            } 
                        }
                    }
                } 
                
                if (breakLoopEM == false || breakLoopType == 2) {
                    HOST_PRINT(commAll, "\n================================================");
                    HOST_PRINT(commAll, "\n========== Start step length search " << equationTypeEM << " ======\n");
                    SLsearchEM.run(commAll, *solverEM, *derivativesEM, receiversEM, sourceSettingsEM, receiversTrueEM, *modelEM, distEM, configEM, modelCoordinatesEM, *gradientEM, steplengthInitEM, dataMisfitEM->getMisfitIt(workflowEM.iteration), workflowEM, freqFilterEM, sourceEstEM, sourceSignalTaperEM, uniqueShotNosRandEM, dataMisfitEM->getMisfitTypeShots());

                    HOST_PRINT(commAll, "================= Update Model " << equationTypeEM << " ============\n\n");
                    /* Apply model update */
                    *gradientEM *= SLsearchEM.getSteplength();
                    *modelEM -= *gradientEM;
                    
                    if (configEM.get<bool>("useModelThresholds"))
                        modelEM->applyThresholds(configEM);

                    if (modelEM->getParameterisation() == 2 || modelEM->getParameterisation() == 1) {
                        HOST_PRINT(commAll, "\n======= calcWaveModulusFromPetrophysics " << equationTypeEM << " =======\n");  
                        modelEM->calcWaveModulusFromPetrophysics();  
                        if (configEM.get<bool>("useModelThresholds"))
                            modelEM->applyThresholds(configEM); 
                    } else if (inversionTypeEM == 3) {
                        HOST_PRINT(commAll, "\n======= calcPetrophysicsFromWaveModulus " << equationTypeEM << " =======\n");  
                        modelEM->calcPetrophysicsFromWaveModulus();
                        if (configEM.get<bool>("useModelThresholds"))
                            modelEM->applyThresholds(configEM); 
                        if (exchangeStrategy != 5 && exchangeStrategy != 6) {
                            // case 0,1,2,3,4: self-constraint of the petrophysical relationship    
                            HOST_PRINT(commAll, "\n======= calcWaveModulusFromPetrophysics " << equationTypeEM << " =======\n");  
                            modelEM->calcWaveModulusFromPetrophysics(); 
                            if (configEM.get<bool>("useModelThresholds"))
                                modelEM->applyThresholds(configEM); 
                        }                    
                    }
                    if (commInterShot->getRank() == 0) {
                        /* only shot Domain 0 writes output */
                        modelEM->write((configEM.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".It_" + std::to_string(workflowEM.iteration + 1)), configEM.get<IndexType>("FileFormat"));
                    }
                    
                    steplengthInitEM *= 0.98;

                    end_t = common::Walltime::get();
                    HOST_PRINT(commAll, "\nFinished iteration " << workflowEM.iteration + 1 << " in " << end_t - start_t << " sec.\n\n");
                }
                
            }  // End of once EM gradient, dataMisfit calculation and model update

            if (inversionTypeEM == 3 && (exchangeStrategy == 1 || exchangeStrategy == 2 || (workflowEM.iteration == maxiterations - 1 && (exchangeStrategy == 3 || exchangeStrategy == 5)))) {
                HOST_PRINT(commAll, "\n=================================================");
                HOST_PRINT(commAll, "\n========= Joint petrophysical inversion =========");
                HOST_PRINT(commAll, "\n=============  From EM to Seismic  ==============");
                HOST_PRINT(commAll, "\n=================================================\n");                
                modelTaper2DJoint.exchangePetrophysics(*modelEM, *model, config);
            }
            
            /* -------------------------------------------------------------------- */
            /* One extra forward modelling to ensure complete and consistent output */
            /* -------------------------------------------------------------------- */
            if (inversionType != 0 && (breakLoop == false || breakLoopType == 2) && workflow.iteration == maxiterations - 1) {
                HOST_PRINT(commAll, "\n================ Maximum number of iterations reached " << equationType << " ================\n");
                HOST_PRINT(commAll, "== Do one more forward modelling to calculate misfit and save seismograms ==\n\n");
            
                if (!useStreamConfig) {                
                    model->prepareForModelling(modelCoordinates, ctx, dist, commShot);
                    solver->prepareForModelling(*model, config.get<ValueType>("DT"));
                }
                                                
                for (IndexType shotInd = shotDist->lb(); shotInd < shotDist->ub(); shotInd++) {
                    if (config.getAndCatch("useRandomSource", 0) == 0) {  
                        shotNumber = uniqueShotNos[shotInd];
                        shotIndTrue = shotInd;
                    } else {
                        shotNumber = uniqueShotNosRand[shotInd];
                        Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNos, shotNumber);
                    }

                    if (useStreamConfig) {
                        HOST_PRINT(commShot, "Switch to model shot: " << shotIndTrue + 1 << " of " << numshots << "\n");
                        model->getModelPerShot(*modelPerShot, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotIndTrue)); 
                        modelPerShot->prepareForModelling(modelCoordinates, ctx, dist, commShot); 
                        solver->prepareForModelling(*modelPerShot, config.get<ValueType>("DT"));
                    }

                    /* Read field data (or pseudo-observed data, respectively) */
                    if (config.get<IndexType>("useReceiversPerShot") == 1) {
                        receivers.init(config, modelCoordinates, ctx, dist, shotNumber);
                        receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber);
                    } else if (config.get<IndexType>("useReceiversPerShot") == 2) {
                        receivers.init(config, modelCoordinates, ctx, dist, shotNumber, numshots);
                        receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber, numshots);
                    }
                    receiversTrue.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".shot_" + std::to_string(shotNumber), 1);

                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshots << ": Additional forward run with " << tStepEnd << " time steps\n");

                    wavefields->resetWavefields();

                    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
                    Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
                    sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);

                    if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0) {
                        sources.getSeismogramHandler().filter(freqFilter);
                        receiversTrue.getSeismogramHandler().filter(freqFilter);
                    }

                    if (config.get<bool>("useSourceSignalInversion")) {
                        sourceEst.applyFilter(sources, shotIndTrue);
                        if (config.get<bool>("useSourceSignalTaper"))
                            sourceSignalTaper.apply(sources.getSeismogramHandler());
                    }

                    if (config.get<IndexType>("useSeismogramTaper") == 2 || config.get<IndexType>("useSeismogramTaper") == 3) {                                                   
                        seismogramTaper2D.apply(receiversTrue.getSeismogramHandler()); 
                    }
                    seismogramTaper1D.apply(receiversTrue.getSeismogramHandler());

                    start_t_shot = common::Walltime::get();

                    if (!useStreamConfig) {
                        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                            solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);
                        }
                    } else {
                        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                            solver->run(receivers, sources, *modelPerShot, *wavefields, *derivatives, tStep);
                        }
                    }

                    // check wavefield and seismogram for NaNs or infinite values
                    if ((commShot->any(!wavefields->isFinite(dist)) || commShot->any(!receivers.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write model and break
                        model->write("model_crash", config.get<IndexType>("FileFormat"));
                        COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield, output model as model_crash.FILE_EXTENSION!");
                    }

                    if (config.get<IndexType>("useSeismogramTaper") == 2 || config.get<IndexType>("useSeismogramTaper") == 3) {                                                   
                        seismogramTaper2D.apply(receivers.getSeismogramHandler()); 
                    }
                    seismogramTaper1D.apply(receivers.getSeismogramHandler());
                        
                    receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
                
                    /* Normalize observed and synthetic data */
                    if (config.get<bool>("normalizeTraces")){
                        receivers.getSeismogramHandler().normalize();
                        receiversTrue.getSeismogramHandler().normalize();
                    }

                    /* Calculate misfit of one shot */
                    misfitPerIt.setValue(shotIndTrue, dataMisfit->calc(receivers, receiversTrue, shotIndTrue));

                    end_t_shot = common::Walltime::get();
                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshots << ": Finished additional forward run\n");

                } //end of loop over Seismic shots
                
                commInterShot->sumArray(misfitPerIt.getLocalValues());
                dataMisfit->addToStorage(misfitPerIt);

                SLsearch.appendToLogFile(commAll, workflow.workflowStage + 1, workflow.iteration + 1, logFilename, dataMisfit->getMisfitSum(workflow.iteration + 1));

                if (workflow.workflowStage != workflow.maxStage - 1) {
                    HOST_PRINT(commAll, "\nChange workflow stage\n");
                    workflow.changeStage(config, *dataMisfit, steplengthInit);
                }
                
            } // end extra Seismic forward modelling

            /* -------------------------------------------------------------------- */
            /* One extra forward modelling to ensure complete and consistent output */
            /* -------------------------------------------------------------------- */
            if (inversionTypeEM != 0 && (breakLoopEM == false || breakLoopType == 2) && workflowEM.iteration == maxiterations - 1) {
                HOST_PRINT(commAll, "\n================ Maximum number of iterations reached " << equationTypeEM << " ================\n");
                HOST_PRINT(commAll, "== Do one more forward modelling to calculate misfit and save seismograms ==\n\n");
            
                if (!useStreamConfigEM) {                
                    modelEM->prepareForModelling(modelCoordinatesEM, ctx, distEM, commShot);
                    solverEM->prepareForModelling(*modelEM, configEM.get<ValueType>("DT"));
                }
                                                
                for (IndexType shotInd = shotDistEM->lb(); shotInd < shotDistEM->ub(); shotInd++) {
                    if (configEM.getAndCatch("useRandomSource", 0) == 0) {  
                        shotNumber = uniqueShotNosEM[shotInd];
                        shotIndTrue = shotInd;
                    } else {
                        shotNumber = uniqueShotNosRandEM[shotInd];
                        Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNosEM, shotNumber);
                    }

                    if (useStreamConfigEM) {
                        HOST_PRINT(commShot, "Switch to model shot: " << shotIndTrue + 1 << " of " << numshotsEM << "\n");
                        modelEM->getModelPerShot(*modelPerShotEM, modelCoordinatesEM, modelCoordinatesBigEM, cutCoordinatesEM.at(shotIndTrue)); 
                        modelPerShotEM->prepareForModelling(modelCoordinatesEM, ctx, distEM, commShot); 
                        solverEM->prepareForModelling(*modelPerShotEM, configEM.get<ValueType>("DT"));
                    }

                    /* Read field data (or pseudo-observed data, respectively) */
                    if (configEM.get<IndexType>("useReceiversPerShot") == 1) {
                        receiversEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber);
                        receiversTrueEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber);
                    } else if (configEM.get<IndexType>("useReceiversPerShot") == 2) {
                        receiversEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber, numshotsEM);
                        receiversTrueEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber, numshotsEM);
                    }
                    receiversTrueEM.getSeismogramHandler().read(configEM.get<IndexType>("SeismogramFormat"), configEM.get<std::string>("fieldSeisName") + ".shot_" + std::to_string(shotNumber), 1);

                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshotsEM << ": Additional forward run with " << tStepEnd << " time steps\n");

                    wavefieldsEM->resetWavefields();

                    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
                    Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettingsEM, shotNumber);
                    sourcesEM.init(sourceSettingsShot, configEM, modelCoordinatesEM, ctx, distEM);

                    if (workflowEM.getLowerCornerFreq() != 0.0 || workflowEM.getUpperCornerFreq() != 0.0) {
                        sourcesEM.getSeismogramHandler().filter(freqFilterEM);
                        receiversTrueEM.getSeismogramHandler().filter(freqFilterEM);
                    }

                    if (configEM.get<bool>("useSourceSignalInversion")) {
                        sourceEstEM.applyFilter(sourcesEM, shotIndTrue);
                        if (configEM.get<bool>("useSourceSignalTaper"))
                            sourceSignalTaperEM.apply(sourcesEM.getSeismogramHandler());
                    }

                    if (configEM.get<IndexType>("useSeismogramTaper") == 2 || configEM.get<IndexType>("useSeismogramTaper") == 3) {                                                   
                        seismogramTaper2DEM.apply(receiversTrueEM.getSeismogramHandler()); 
                    }
                    seismogramTaper1DEM.apply(receiversTrueEM.getSeismogramHandler());

                    start_t_shot = common::Walltime::get();

                    if (!useStreamConfigEM) {
                        for (IndexType tStep = 0; tStep < tStepEndEM; tStep++) {
                            solverEM->run(receiversEM, sourcesEM, *modelEM, *wavefieldsEM, *derivativesEM, tStep);
                        }
                    } else {
                        for (IndexType tStep = 0; tStep < tStepEndEM; tStep++) {
                            solverEM->run(receiversEM, sourcesEM, *modelPerShotEM, *wavefieldsEM, *derivativesEM, tStep);
                        }
                    }

                    // check wavefield and seismogram for NaNs or infinite values
                    if ((commShot->any(!wavefieldsEM->isFinite(distEM)) || commShot->any(!receiversEM.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write modelEM and break
                        modelEM->write("model_crash", configEM.get<IndexType>("FileFormat"));
                        COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield, output model as model_crash.FILE_EXTENSION!");
                    }

                    if (configEM.get<IndexType>("useSeismogramTaper") == 2 || configEM.get<IndexType>("useSeismogramTaper") == 3) {                                                   
                        seismogramTaper2DEM.apply(receiversEM.getSeismogramHandler()); 
                    }
                    seismogramTaper1DEM.apply(receiversEM.getSeismogramHandler());
                        
                    receiversEM.getSeismogramHandler().write(configEM.get<IndexType>("SeismogramFormat"), configEM.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".It_" + std::to_string(workflowEM.iteration + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinatesEM);
                
                    /* Normalize observed and synthetic data */
                    if (configEM.get<bool>("normalizeTraces")){
                        receiversEM.getSeismogramHandler().normalize();
                        receiversTrueEM.getSeismogramHandler().normalize();
                    }

                    /* Calculate misfit of one shot */
                    misfitPerItEM.setValue(shotIndTrue, dataMisfitEM->calc(receiversEM, receiversTrueEM, shotIndTrue));

                    end_t_shot = common::Walltime::get();
                    HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshotsEM << ": Finished additional forward run\n");

                } //end of loop over EM shots
                commInterShot->sumArray(misfitPerItEM.getLocalValues());
                dataMisfitEM->addToStorage(misfitPerItEM);

                SLsearchEM.appendToLogFile(commAll, workflowEM.workflowStage + 1, workflowEM.iteration + 1, logFilenameEM, dataMisfitEM->getMisfitSum(workflowEM.iteration + 1));

                if (workflowEM.workflowStage != workflowEM.maxStage - 1) {
                    HOST_PRINT(commAll, "\nChange workflowEM stage\n");
                    workflowEM.changeStage(configEM, *dataMisfitEM, steplengthInitEM);
                }
                
            } // end extra EM forward modelling

        } // end of loop over iterations 
        
    } // end of loop over workflow stages 
    
    if (inversionType == 3 && (exchangeStrategy == 4 || exchangeStrategy == 6)) {   
        HOST_PRINT(commAll, "\n=================================================");
        HOST_PRINT(commAll, "\n========= Joint petrophysical inversion =========");
        HOST_PRINT(commAll, "\n=============  From Seismic to EM  ==============");
        HOST_PRINT(commAll, "\n=================================================\n");
        modelTaper2DJoint.exchangePetrophysics(*model, *modelEM, configEM); 
        if (commInterShot->getRank() == 0) {
            /* only shot Domain 0 writes output */
            if (!useStreamConfigEM) {
                modelEM->write(configEM.get<std::string>("ModelFilename"), configEM.get<IndexType>("FileFormat"));
            } else {
                modelEM->write(configBigEM.get<std::string>("ModelFilename"), configEM.get<IndexType>("FileFormat"));
            }
        }
    }
    
    globalEnd_t = common::Walltime::get();
    HOST_PRINT(commAll, "\nTotal runtime of WAVE-Inversion: " << globalEnd_t - globalStart_t << " sec.\nWAVE-Inversion finished!\n\n");
    return 0;
}
