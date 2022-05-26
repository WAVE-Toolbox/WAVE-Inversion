
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

#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolver/ForwardSolverFactory.hpp>
#include <ForwardSolver/SourceReceiverImpl/SourceReceiverImpl.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>

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
#include "Common/InversionSingle.hpp"

using namespace scai;
using namespace KITGPI;

bool verbose; // global variable definition

int main(int argc, char *argv[])
{
    double globalStart_t, globalEnd_t; /* For timing */
    globalStart_t = common::Walltime::get();
    IndexType seedtime = (int)time(0);
    
    /* --------------------------------------- */
    /* Read configuration from file            */
    /* --------------------------------------- */
    Configuration::Configuration config;
    Configuration::Configuration configEM; 
    IndexType inversionType = 0;
    IndexType inversionTypeEM = 0;
    if (argc == 2) {
        config.readFromFile(argv[1], true);
        configEM.readFromFile(argv[1], true);
        inversionType = config.getAndCatch("inversionType", 1);
    } else if (argc == 3) {
        config.readFromFile(argv[1], true);
        configEM.readFromFile(argv[2], true);
        inversionType = config.getAndCatch("inversionType", 1);
        inversionTypeEM = configEM.getAndCatch("inversionType", 1);
        if (inversionType == 0)
            config.readFromFile(argv[2], true);
        if (inversionTypeEM == 0)
            configEM.readFromFile(argv[1], true);
    } else {
        std::cout << "\n\nNo configuration file given!\n\n" << std::endl;
        return (2);
    }
    verbose = config.get<bool>("verbose");
    
    // Initialization of inversion
    /* inversionType: 0 = no inversion, 1 = single inversion, 2 = joint structural inversion, 3 = joint petrophysical inversion, 4 = joint inversion of the seismic waves or EM waves */    
    InversionSingle<ValueType> inversionSingle(config, inversionType);
    InversionSingle<ValueType> inversionSingleEM(configEM, inversionTypeEM);
    
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);   
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);  
    bool isSeismic = Common::checkEquationType<ValueType>(equationType);
    
    std::string dimensionEM = configEM.get<std::string>("dimension");
    std::string equationTypeEM = configEM.get<std::string>("equationType"); 
    std::transform(dimensionEM.begin(), dimensionEM.end(), dimensionEM.begin(), ::tolower);   
    std::transform(equationTypeEM.begin(), equationTypeEM.end(), equationTypeEM.begin(), ::tolower);  
    bool isSeismicEM = Common::checkEquationType<ValueType>(equationTypeEM);      
        
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
    
    /* exchangeStrategy: 0 = self-restraint, 1 = self-restraint + mutual restraint (single parameter), 2 = self-restraint + mutual restraint (all parameters), 3 = self-restraint + mutual restraint + staged exchange, 4 = self-restraint + mutual restraint + sequential exchange, 5 = mutual restraint + staged exchange, 6 = mutual restraint + sequential exchange. */
    IndexType exchangeStrategy = config.get<IndexType>("exchangeStrategy");
    IndexType exchangeStrategyEM = configEM.get<IndexType>("exchangeStrategy");
    SCAI_ASSERT_ERROR(exchangeStrategy == exchangeStrategyEM, "exchangeStrategy != exchangeStrategyEM"); // check whether exchangeStrategy has been applied successfully.
        
    /*breakLoop: 0 = break iteration loop of which the abort criterion is true while do not influence another one. 0 is the same as individual inversion. 1 = break iteration loop if one of the two abort criterion is true. 2 = break iteration loop if both the two abort criterion is true. */
    IndexType breakLoopType = config.get<IndexType>("breakLoopType"); 
    IndexType breakLoopTypeEM = configEM.get<IndexType>("breakLoopType");
    SCAI_ASSERT_ERROR(breakLoopType == breakLoopTypeEM, "breakLoopType != breakLoopTypeEM"); // check whether breakLoopType has been applied successfully.
    if (inversionType == 1 || inversionTypeEM == 1 || exchangeStrategy > 2) {
        SCAI_ASSERT_ERROR(breakLoopType == 0, "breakLoopType != 0"); // check whether breakLoopType has been applied successfully.
    }
    if (inversionType == 4 && inversionTypeEM == 4) {
        SCAI_ASSERT_ERROR(isSeismic == isSeismicEM, "isSeismic != isSeismicEM");
    }

    std::string misfitType = config.get<std::string>("misfitType");
    std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::tolower);
    std::string logFilename = config.get<std::string>("logFilename");
    IndexType gradientKernel = config.getAndCatch("gradientKernel", 0);
    
    std::string misfitTypeEM = configEM.get<std::string>("misfitType");
    std::transform(misfitTypeEM.begin(), misfitTypeEM.end(), misfitTypeEM.begin(), ::tolower);
    std::string logFilenameEM = configEM.get<std::string>("logFilename");
    IndexType gradientKernelEM = configEM.getAndCatch("gradientKernel", 0);
    
    IndexType maxiterations = config.get<IndexType>("maxIterations");
    if (inversionType == 0 && inversionTypeEM != 0) {
        maxiterations = configEM.get<IndexType>("maxIterations");
    } 
    
    /* inter node communicator */
    dmemo::CommunicatorPtr commAll = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(commAll->getNodeRank());

    std::string settingsFilename; // filename for processor specific settings
    if (common::Settings::getEnvironment(settingsFilename, "SCAI_SETTINGS")) {
        // each processor reads line of settings file that matches its node name and node rank
        common::Settings::readSettingsFile(settingsFilename.c_str(), commAll->getNodeName(), commAll->getNodeRank());
    }    
    
    inversionSingle.printConfig(commAll, config, inversionType, 1);
    inversionSingleEM.printConfig(commAll, configEM, inversionTypeEM, 2);
    
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    // distribute the grid onto available processors
    dmemo::DistributionPtr dist = nullptr;
    dmemo::DistributionPtr distEM = nullptr;
    dmemo::DistributionPtr distBig = nullptr;
    dmemo::DistributionPtr distBigEM = nullptr;
    
    /* --------------------------------------- */
    /* coordinate mapping (3D<->1D)            */
    /* --------------------------------------- */
    Acquisition::Coordinates<ValueType> modelCoordinates;
    Acquisition::Coordinates<ValueType> modelCoordinatesBig;
    Acquisition::Coordinates<ValueType> modelCoordinatesEM;
    Acquisition::Coordinates<ValueType> modelCoordinatesBigEM;
    
    /* --------------------------------------- */
    /* Factories                               */
    /* --------------------------------------- */    
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));    
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr modelEM(Modelparameter::Factory<ValueType>::Create(equationTypeEM));
    
    /* --------------------------------------- */
    /* Calculate derivative matrices           */
    /* --------------------------------------- */
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivativesInversion(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivativesInversionEM(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimensionEM));
            
    /* --------------------------------------- */
    /* Misfit                                  */
    /* --------------------------------------- */
    Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(Misfit::Factory<ValueType>::Create(misfitType));    
    Misfit::Misfit<ValueType>::MisfitPtr dataMisfitEM(Misfit::Factory<ValueType>::Create(misfitTypeEM));
    
    /* --------------------------------------- */
    /* Workflow                                */
    /* --------------------------------------- */
    Workflow::Workflow<ValueType> workflow;
    Workflow::Workflow<ValueType> workflowEM;
    workflow.init(config);
    workflowEM.init(configEM);

    /* --------------------------------------- */
    /* Abort criterion                         */
    /* --------------------------------------- */
    bool breakLoop = false;
    bool breakLoopEM = false;
    StepLengthSearch<ValueType> SLsearch;
    StepLengthSearch<ValueType> SLsearchEM;

    /* --------------------------------------- */
    /* Gradients                               */
    /* --------------------------------------- */
    typename Gradient::Gradient<ValueType>::GradientPtr crossGradientDerivative(Gradient::Factory<ValueType>::Create(equationType));    
    typename Gradient::Gradient<ValueType>::GradientPtr  crossGradientDerivativeEM(Gradient::Factory<ValueType>::Create(equationTypeEM)); 
            
    /* --------------------------------------- */
    /* Gradient taper                          */
    /* --------------------------------------- */
    Taper::Taper2D<ValueType> modelTaper2DJoint;    
    
    inversionSingle.init(commAll, config, inversionType, 1, modelCoordinates, modelCoordinatesBig, ctx, dist, distBig, derivativesInversion, maxiterations, model, workflow, crossGradientDerivative, SLsearch);
    inversionSingleEM.init(commAll, configEM, inversionTypeEM, 2, modelCoordinatesEM, modelCoordinatesBigEM, ctx, distEM, distBigEM, derivativesInversionEM, maxiterations, modelEM, workflowEM, crossGradientDerivativeEM, SLsearchEM);
    
    inversionSingle.estimateMemory(commAll, config, inversionType, 1, dist, modelCoordinates, model);
    inversionSingleEM.estimateMemory(commAll, configEM, inversionTypeEM, 2, distEM, modelCoordinatesEM, modelEM);
         
    // to ensure the self-constraint of individual FWI, e.g., structural constraint of vs on density in SH FWI.
    if (inversionType != 0 && inversionTypeEM == 0) {
        *modelEM = *model;
        *derivativesInversionEM = *derivativesInversion;
        distEM = dist;
        distBigEM = distBig;
        modelCoordinatesBigEM = modelCoordinatesBig;
    } else if (inversionType == 0 && inversionTypeEM != 0) {
        *model = *modelEM;
        *derivativesInversion = *derivativesInversionEM;
        dist = distEM;
        distBig = distBigEM;
        modelCoordinatesBig = modelCoordinatesBigEM;
    }
    if (inversionType > 1 || inversionTypeEM > 1) {
        if (!useStreamConfig && !useStreamConfigEM) {
            modelTaper2DJoint.initTransformMatrix(dist, distEM, ctx);  
            modelTaper2DJoint.calcTransformMatrix1to2(modelCoordinates, modelCoordinatesEM);
            modelTaper2DJoint.calcTransformMatrix2to1(modelCoordinates, modelCoordinatesEM); 
        } else if (useStreamConfig && !useStreamConfigEM) {
            modelTaper2DJoint.initTransformMatrix(distBig, distEM, ctx);  
            modelTaper2DJoint.calcTransformMatrix1to2(modelCoordinatesBig, modelCoordinatesEM);
            modelTaper2DJoint.calcTransformMatrix2to1(modelCoordinatesBig, modelCoordinatesEM); 
        } else if (!useStreamConfig && useStreamConfigEM) {
            modelTaper2DJoint.initTransformMatrix(dist, distBigEM, ctx);  
            modelTaper2DJoint.calcTransformMatrix1to2(modelCoordinates, modelCoordinatesBigEM);
            modelTaper2DJoint.calcTransformMatrix2to1(modelCoordinates, modelCoordinatesBigEM); 
        } else if (useStreamConfig && useStreamConfigEM) {
            modelTaper2DJoint.initTransformMatrix(distBig, distBigEM, ctx);  
            modelTaper2DJoint.calcTransformMatrix1to2(modelCoordinatesBig, modelCoordinatesBigEM);
            modelTaper2DJoint.calcTransformMatrix2to1(modelCoordinatesBig, modelCoordinatesBigEM); 
        }
    }
    
    globalEnd_t = common::Walltime::get();
    HOST_PRINT(commAll, "\nFinished all initialization in " << globalEnd_t - globalStart_t << " sec.\n");
    
    /* --------------------------------------- */
    /*       Loop over workflow stages         */
    /* --------------------------------------- */
    IndexType stageCount = 0;
    for (workflow.workflowStage = 0; workflow.workflowStage < workflow.maxStage; workflow.workflowStage++) {
        workflowEM.workflowStage = workflow.workflowStage;
        bool breakLoopLast = breakLoop;
        bool breakLoopLastEM = breakLoopEM;
        IndexType useRTM = 0; 
        IndexType useRTMEM = 0; 
        
        inversionSingle.initStage(commAll, config, inversionType, 1, ctx, workflow, dataMisfit, breakLoop, dist, breakLoopEM);        
        inversionSingleEM.initStage(commAll, configEM, inversionTypeEM, 2, ctx, workflowEM, dataMisfitEM, breakLoopEM, distEM, breakLoop);
        
        /* --------------------------------------- */
        /*        Loop over iterations             */
        /* --------------------------------------- */ 
        for (workflow.iteration = 0; workflow.iteration < maxiterations; workflow.iteration++) {
            workflowEM.iteration = workflow.iteration;
            /* --------------------------------------- */
            /*        Start the first inversion        */
            /* --------------------------------------- */   
            inversionSingle.calcGradient(commAll, dist, model, config, modelCoordinates, modelCoordinatesBig, workflow, dataMisfit, crossGradientDerivative, SLsearch, modelTaper2DJoint, maxiterations, useRTM, breakLoop, ctx, seedtime, inversionType, 1, breakLoopEM, modelEM, configEM, modelCoordinatesEM, workflowEM, dataMisfitEM, crossGradientDerivativeEM, derivativesInversionEM);
            
            if (useRTM == 1) {
                HOST_PRINT(commAll, "\nFinish inverse time migration after stage " << workflow.workflowStage + 1 << " in " << equationType << " 1 \n");
                break;
            }
            
            // We set a new break condition so that two inversions can change stage simultaneously in joint inversion
            if (breakLoop == true) {                 
                if (inversionType > 2 && (exchangeStrategy == 1 || exchangeStrategy == 2 || exchangeStrategy == 3 || exchangeStrategy == 5)) { 
                    if (inversionType == 3) {
                        modelTaper2DJoint.exchangePetrophysics(commAll, *model, config, *modelEM, configEM, 2); 
                    } else if (inversionType == 4) {
                        modelTaper2DJoint.exchangeModelparameters(commAll, *model, config, *modelEM, configEM, 2); 
                    }            
                }       
                if (breakLoopType == 0 && breakLoopEM == true) {
                    if (exchangeStrategy == 3 || exchangeStrategy == 5) {
                        breakLoopEM = false;
                        workflow.iteration = 0;
                    } else if (gradientKernel == 4) {
                        useRTM = 1;
                    } else {
                        break;
                    }                       
                }           
            } 

            inversionSingle.updateModel(commAll, dist, model, config, modelCoordinates, modelCoordinatesBig, workflow, dataMisfit, crossGradientDerivative, SLsearch, modelTaper2DJoint, maxiterations, useRTM, breakLoop, ctx, seedtime, inversionType, 1, breakLoopEM, modelEM, configEM, modelCoordinatesEM, workflowEM, dataMisfitEM, crossGradientDerivativeEM, derivativesInversionEM);
            
            if (inversionType > 2 && (breakLoop == false || breakLoopType == 2) && (exchangeStrategy == 1 || exchangeStrategy == 2 || (workflow.iteration == maxiterations - 1 && (exchangeStrategy == 3 || exchangeStrategy == 5)))) {
                if (inversionType == 3) {
                    modelTaper2DJoint.exchangePetrophysics(commAll, *model, config, *modelEM, configEM, 2); 
                } else if (inversionType == 4) {
                    modelTaper2DJoint.exchangeModelparameters(commAll, *model, config, *modelEM, configEM, 2); 
                }            
            }            
            
            inversionSingle.runExtra(commAll, dist, model, config, modelCoordinates, modelCoordinatesBig, workflow, dataMisfit, crossGradientDerivative, SLsearch, modelTaper2DJoint, maxiterations, useRTM, breakLoop, ctx, seedtime, inversionType, 1, breakLoopEM, modelEM, configEM, modelCoordinatesEM, workflowEM, dataMisfitEM, crossGradientDerivativeEM, derivativesInversionEM);
            
            /* --------------------------------------- */
            /*        Start the second inversion       */
            /* --------------------------------------- */            
            inversionSingleEM.calcGradient(commAll, distEM, modelEM, configEM, modelCoordinatesEM, modelCoordinatesBigEM, workflowEM, dataMisfitEM, crossGradientDerivativeEM, SLsearchEM, modelTaper2DJoint, maxiterations, useRTMEM, breakLoopEM, ctx, seedtime, inversionTypeEM, 2, breakLoop, model, config, modelCoordinates, workflow, dataMisfit, crossGradientDerivative, derivativesInversion);   
           
            if (useRTMEM == 1) {
                HOST_PRINT(commAll, "\nFinish inverse time migration after stage " << workflowEM.workflowStage + 1 << " in " << equationTypeEM << " 2 \n");
                break;
            }

            // We set a new break condition so that two inversions can change stage simultaneously in joint inversion
            if (breakLoopEM == true) {   
                if (inversionTypeEM > 2 && (exchangeStrategy == 1 || exchangeStrategy == 2 || exchangeStrategy == 3 || exchangeStrategy == 5)) {
                    if (inversionTypeEM == 3) {
                        modelTaper2DJoint.exchangePetrophysics(commAll, *modelEM, configEM, *model, config, 1); 
                    } else if (inversionTypeEM == 4) {
                        modelTaper2DJoint.exchangeModelparameters(commAll, *modelEM, configEM, *model, config, 1); 
                    } 
                }   
                if (breakLoopType == 1) {
                    if (gradientKernelEM == 4) {
                        useRTMEM = 1;
                    } else {
                        break;
                    }       
                } else if (breakLoopType == 2) {
                    if (breakLoop == true) {
                        if (gradientKernelEM == 4) {
                            useRTMEM = 1;
                        } else {
                            break;
                        }       
                    } else {
                        breakLoopEM = false;
                    }
                } else if (breakLoopType == 0) {
                    if (gradientKernelEM == 4) {
                        useRTMEM = 1;
                    }
                    if (breakLoop == true) {
                        if (gradientKernelEM != 4) {
                            break;
                        }    
                    } else {
                        if (inversionType > 2 && (exchangeStrategy == 1 || exchangeStrategy == 2 || exchangeStrategy == 3 || exchangeStrategy == 5)) {
                            if (inversionType == 3) {
                                modelTaper2DJoint.exchangePetrophysics(commAll, *model, config, *modelEM, configEM, 2); 
                            } else if (inversionType == 4) {
                                modelTaper2DJoint.exchangeModelparameters(commAll, *model, config, *modelEM, configEM, 2); 
                            }           
                        } 
                    }
                }
            } 
 
            inversionSingleEM.updateModel(commAll, distEM, modelEM, configEM, modelCoordinatesEM, modelCoordinatesBigEM, workflowEM, dataMisfitEM, crossGradientDerivativeEM, SLsearchEM, modelTaper2DJoint, maxiterations, useRTMEM, breakLoopEM, ctx, seedtime, inversionTypeEM, 2, breakLoop, model, config, modelCoordinates, workflow, dataMisfit, crossGradientDerivative, derivativesInversionEM);    
               
            if (inversionTypeEM > 2 && (breakLoopEM == false || breakLoopType == 2) && (exchangeStrategy == 1 || exchangeStrategy == 2 || (workflowEM.iteration == maxiterations - 1 && (exchangeStrategy == 3 || exchangeStrategy == 5)))) {
                if (inversionTypeEM == 3) {
                    modelTaper2DJoint.exchangePetrophysics(commAll, *modelEM, configEM, *model, config, 1); 
                } else if (inversionTypeEM == 4) {
                    modelTaper2DJoint.exchangeModelparameters(commAll, *modelEM, configEM, *model, config, 1); 
                }           
            }

            inversionSingleEM.runExtra(commAll, distEM, modelEM, configEM, modelCoordinatesEM, modelCoordinatesBigEM, workflowEM, dataMisfitEM, crossGradientDerivativeEM, SLsearchEM, modelTaper2DJoint, maxiterations, useRTMEM, breakLoopEM, ctx, seedtime, inversionTypeEM, 2, breakLoop, model, config, modelCoordinates, workflow, dataMisfit, crossGradientDerivative, derivativesInversionEM);  
        } // end of loop over iterations 
    
        if (workflow.workflowStage < workflow.maxStage - 1) {
            if (breakLoop == true && breakLoopEM == true) {
                if (exchangeStrategy == 1 || exchangeStrategy == 2 || exchangeStrategy == 3 || exchangeStrategy == 5) {
                    breakLoop = false;
                    breakLoopEM = false;
                } else if (exchangeStrategy == 4 || exchangeStrategy == 6) {
                    if (breakLoop != breakLoopLast)
                        breakLoop = false;
                    if (breakLoopEM != breakLoopLastEM)
                        breakLoopEM = false;
                }
            } 
            if (inversionType == 1 || inversionTypeEM == 1 || exchangeStrategy == 0) {
                if (breakLoop != breakLoopLast)
                    breakLoop = false;
                if (breakLoopEM != breakLoopLastEM)
                    breakLoopEM = false;
            }
        }
        
        if (workflow.workflowStage == workflow.maxStage - 1 && inversionType > 1 && (exchangeStrategy == 4 || exchangeStrategy == 6)) { 
            stageCount++;
            if (inversionType == 3) {
                modelTaper2DJoint.exchangePetrophysics(commAll, *model, config, *modelEM, configEM, 2); 
            } else if (inversionType == 4) {
                if (stageCount % 2 == 1) {
                    modelTaper2DJoint.exchangeModelparameters(commAll, *model, config, *modelEM, configEM, 2); 
                    HOST_PRINT(commAll, "\nChange workflow stage from " << equationType << " 1 to " << equationTypeEM << " 2\n");
                } else {
                    modelTaper2DJoint.exchangeModelparameters(commAll, *modelEM, configEM, *model, config, 1); 
                    
                    SLsearch.initLogFile(commAll, logFilename, misfitType, config.getAndCatch("steplengthType", 2), workflow.getInvertForParameters().size(), config.getAndCatch("saveCrossGradientMisfit", 0));
                    HOST_PRINT(commAll, "\nChange workflow stage from " << equationTypeEM << " 2 to " << equationType << " 1\n");
                }
            } 
            if ((inversionType == 3 && stageCount < 2) || (inversionType == 4 && stageCount < 3)) {
                workflow.workflowStage = -1;
                workflow.init(config);
                workflowEM.init(configEM);
            }
            if (stageCount % 2 == 1) {
                breakLoop = true;
                breakLoopEM = false;
            } else if (breakLoopEM != breakLoopLastEM) {
                breakLoop = false;
                breakLoopEM = true;
            }
        }    
    } // end of loop over workflow stages 
    
    globalEnd_t = common::Walltime::get();
    if (inversionType != 0)
        SLsearch.appendRunTimeToLogFile(commAll, logFilename, globalEnd_t - globalStart_t);
    if (inversionTypeEM != 0)
        SLsearchEM.appendRunTimeToLogFile(commAll, logFilenameEM, globalEnd_t - globalStart_t);
    HOST_PRINT(commAll, "\nTotal runtime of WAVE-Inversion: " << globalEnd_t - globalStart_t << " sec.\nWAVE-Inversion finished!\n\n");
    return 0;
}
