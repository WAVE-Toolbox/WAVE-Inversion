#include <scai/lama.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/dmemo/CommunicatorStack.hpp>

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include <Configuration/Configuration.hpp>
//#include "/home/tmetz/projects/FDSimulation_LAMA//src/Configuration/Configuration.hpp"
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>

#include <ForwardSolver/ForwardSolver.hpp>

#include <CheckParameter/CheckParameter.hpp>
#include <Filter/Filter.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolver/ForwardSolverFactory.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>

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
#include "Taper/TaperFactory.hpp"
#include "Workflow/Workflow.hpp"

#include <Common/HostPrint.hpp>
#include <Partitioning/PartitioningCubes.hpp>

using namespace scai;
using namespace KITGPI;

bool verbose; // global variable definition

int main(int argc, char *argv[])
{
    typedef double ValueType;
    double start_t, end_t, start_t_shot, end_t_shot; /* For timing */

    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    /* --------------------------------------- */
    /* Read configuration from file            */
    /* --------------------------------------- */
    Configuration::Configuration config(argv[1]);

    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");

    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    std::string misfitType = config.get<std::string>("misfitType");
    std::string fieldSeisName(config.get<std::string>("fieldSeisName"));
    std::string gradname(config.get<std::string>("gradientFilename"));
    std::string logFilename = config.get<std::string>("logFilename");
    ValueType steplengthInit = config.get<ValueType>("steplengthInit");
    IndexType maxiterations = config.get<IndexType>("maxIterations");
    std::string optimizationType = config.get<std::string>("optimizationType");

    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */
    /* inter node communicator */
    dmemo::CommunicatorPtr commAll = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(commAll->getNodeRank());
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    // inter node distribution
    // define the grid topology by sizes NX, NY, and NZ from configuration
    // Attention: LAMA uses row-major indexing while SOFI-3D uses column-major, so switch dimensions, x-dimension has stride 1

    IndexType npS = config.get<IndexType>("ProcNS");
    IndexType npX = config.get<IndexType>("ProcNX");
    IndexType npY = config.get<IndexType>("ProcNY");
    IndexType npZ = config.get<IndexType>("ProcNZ");

    // following lines should be part of checkParameter.tpp
    if (commAll->getSize() != npS * npX * npY * npZ) {
        HOST_PRINT(commAll, "\n Error: Number of MPI processes (" << commAll->getSize() << ") doesn't match the number of processes specified in " << argv[1] << ": ProcNS * ProcNX * ProcNY * ProcNZ = " << npS * npX * npY * npZ << "\n")
        return (2);
    }

    // Build subsets of processors for the shots

    common::Grid2D procAllGrid(npS, npX * npY * npZ);
    IndexType procAllGridRank[2];

    procAllGrid.gridPos(procAllGridRank, commAll->getRank());

    // communicator for set of processors that solve one shot

    dmemo::CommunicatorPtr commShot = commAll->split(procAllGridRank[0]);

    // this communicator is used for reducing the solutions of problems

    dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    SCAI_DMEMO_TASK(commShot)

    IndexType nx = config.get<IndexType>("NX");
    IndexType ny = config.get<IndexType>("NY");
    IndexType nz = config.get<IndexType>("NZ");

    common::Grid3D grid(nz, ny, nx);
    common::Grid3D procGrid(npZ, npY, npX);
    // distribute the grid onto available processors, topology can be set by environment variable
    dmemo::DistributionPtr dist(new dmemo::GridDistribution(grid, commShot, procGrid));

    verbose = config.get<bool>("verbose");
    HOST_PRINT(commAll, "\nIFOS" << dimension << " " << equationType << " - LAMA Version\n\n");
    if (commAll->getRank() == MASTERGPI) {
        config.print();
    }

    /* --------------------------------------- */
    /* Calculate derivative matrices           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    derivatives->init(dist, ctx, config, commShot);
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing matrices in " << end_t - start_t << " sec.\n\n");

    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    Acquisition::Sources<ValueType> sources(config, ctx, dist);
    CheckParameter::checkSources<ValueType>(config, sources, commShot);
    dmemo::BlockDistribution shotDist(sources.getNumShots(), commInterShot);
    Acquisition::Receivers<ValueType> receivers;
    if (!config.get<bool>("useReceiversPerShot"))
        receivers.init(config, ctx, dist);

    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    // load starting model
    model->init(config, ctx, dist);

    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields = Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefields->init(ctx, dist);

    //Temporary Wavefield for the derivative of the forward wavefields
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldsTemp = Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefieldsTemp->init(ctx, dist);

    /* --------------------------------------- */
    /* Wavefield record                        */
    /* --------------------------------------- */
    typedef typename Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
    std::vector<wavefieldPtr> wavefieldrecord;

    for (IndexType i = 0; i < tStepEnd; i++) {
        wavefieldPtr wavefieldsTemp(Wavefields::Factory<ValueType>::Create(dimension, equationType));
        wavefieldsTemp->init(ctx, dist);

        wavefieldrecord.push_back(wavefieldsTemp);
    }

    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */
    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));
    solver->initForwardSolver(config, *derivatives, *wavefields, *model, ctx, config.get<ValueType>("DT"));

    /* --------------------------------------- */
    /* True data                               */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> receiversTrue;
    if (!config.get<bool>("useReceiversPerShot")) {
        receiversTrue.init(config, ctx, dist);
        CheckParameter::checkReceivers<ValueType>(config, receiversTrue, commShot);
    }

    /* --------------------------------------- */
    /* Misfit                                  */
    /* --------------------------------------- */
    Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(Misfit::Factory<ValueType>::Create(misfitType));
    lama::DenseVector<ValueType> misfitPerIt(sources.getNumShots(), 0, ctx);

    /* --------------------------------------- */
    /* Adjoint sources                         */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> adjointSources;
    if (!config.get<bool>("useReceiversPerShot"))
        adjointSources.init(config, ctx, dist);

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
    Taper::Taper<ValueType>::TaperPtr sourceSignalTaper(Taper::Factory<ValueType>::Create("1D"));
    if (config.get<bool>("useSourceSignalInversion")) {
        sourceEst.init(tStepEnd, sources.getCoordinates().getDistributionPtr(), config.get<ValueType>("waterLevel"));
        if (config.get<bool>("useSourceSignalTaper")) {
            sourceSignalTaper->init(std::make_shared<dmemo::NoDistribution>(tStepEnd), ctx, 1);
            if (config.get<IndexType>("sourceSignalTaperStart2") == 0 && config.get<IndexType>("sourceSignalTaperEnd2") == 0)
                sourceSignalTaper->calcCosineTaper(config.get<IndexType>("sourceSignalTaperStart1"), config.get<IndexType>("sourceSignalTaperEnd1"), 0);
            else
                sourceSignalTaper->calcCosineTaper(config.get<IndexType>("sourceSignalTaperStart1"), config.get<IndexType>("sourceSignalTaperEnd1"), config.get<IndexType>("sourceSignalTaperStart2"), config.get<IndexType>("sourceSignalTaperEnd2"), 0);
        }
    }

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

    /* --------------------------------------- */
    /* Gradient calculation                    */
    /* --------------------------------------- */
    GradientCalculation<ValueType> gradientCalculation;

    /* --------------------------------------- */
    /* Gradient taper                          */
    /* --------------------------------------- */
    Preconditioning::SourceReceiverTaper<ValueType> ReceiverTaper;
    if (!config.get<bool>("useReceiversPerShot"))
        ReceiverTaper.init(dist, ctx, receivers, config, config.get<IndexType>("receiverTaperRadius"));
    Taper::Taper<ValueType>::TaperPtr gradientTaper(Taper::Factory<ValueType>::Create("1D"));
    if (config.get<bool>("useGradientTaper")) {
        gradientTaper->init(dist, ctx, 1);
        gradientTaper->readTaper(config.get<std::string>("gradientTaperName") + ".mtx", config.get<IndexType>("PartitionedIn"));
    }

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

        /* --------------------------------------- */
        /*        Loop over iterations             */
        /* --------------------------------------- */

        for (workflow.iteration = 0; workflow.iteration < maxiterations; workflow.iteration++) {

            HOST_PRINT(commAll, "\n=================================================");
            HOST_PRINT(commAll, "\n============ Workflow stage " << workflow.workflowStage + 1 << " of " << workflow.maxStage << " ==============");
            HOST_PRINT(commAll, "\n============     Iteration " << workflow.iteration + 1 << "       ==============");
            HOST_PRINT(commAll, "\n=================================================\n\n");
            start_t = common::Walltime::get();

            /* Update model for fd simulation (averaging, inverse Density ...) */
            model->prepareForModelling(config, ctx, dist, commShot);
            solver->prepareForModelling(*model, config.get<ValueType>("DT"));

            /* --------------------------------------- */
            /*        Loop over shots                  */
            /* --------------------------------------- */

            gradient->resetGradient(); // reset gradient because gradient is a sum of all gradientsPerShot gradients+=gradientPerShot
            misfitPerIt = 0;

            for (IndexType shotNumber = shotDist.lb(); shotNumber < shotDist.ub(); shotNumber++) {

                HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << sources.getNumShots() << ": started\n");

                if (config.get<bool>("useReceiversPerShot")) {
                    receivers.init(config, ctx, dist, shotNumber);
                    receiversTrue.init(config, ctx, dist, shotNumber);
                    adjointSources.init(config, ctx, dist, shotNumber);

                    CheckParameter::checkReceivers<ValueType>(config, receivers, commShot);

                    ReceiverTaper.init(dist, ctx, receivers, config, config.get<IndexType>("receiverTaperRadius"));
                }
                /* Read field data (or pseudo-observed data, respectively) */
                receiversTrue.getSeismogramHandler().read(config, fieldSeisName + ".shot_" + std::to_string(shotNumber), 1);
                if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0)
                    receiversTrue.getSeismogramHandler().filter(freqFilter);

                /* Reset approximated Hessian per shot */
                if (config.get<bool>("useEnergyPreconditioning") == 1)
                    energyPrecond.resetApproxHessian();

                sources.init(config, ctx, dist, shotNumber);
                if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0)
                    sources.getSeismogramHandler().filter(freqFilter);

                /* Source time function inversion */
                if (config.get<bool>("useSourceSignalInversion") == 1) {
                    if (workflow.iteration == 0) {
                        HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << sources.getNumShots() << " : Source Time Function Inversion\n");

                        wavefields->resetWavefields();

                        for (scai::IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                            solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);
                        }

                        solver->resetCPML();

                        if (config.get<bool>("maxOffsetSrcEst") == 1)
                            sourceEst.calcOffsetMutes(sources, receivers, config.get<ValueType>("maxOffsetSrcEst"), nx, ny, nz);

                        sourceEst.estimateSourceSignal(receivers, receiversTrue, shotNumber);
                        sourceEst.applyFilter(sources, shotNumber);
                        if (config.get<bool>("useSourceSignalTaper"))
                            sourceSignalTaper->apply(sources.getSeismogramHandler());

                        if (config.get<bool>("writeInvertedSource") == 1)
                            sources.getSeismogramHandler().write(config, config.get<std::string>("sourceSeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber));
                    } else {
                        sourceEst.applyFilter(sources, shotNumber);
                        sourceSignalTaper->apply(sources.getSeismogramHandler());
                    }
                }

                /* --------------------------------------- */
                /*        Forward modelling                */
                /* --------------------------------------- */

                HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << sources.getNumShots() << ": Start time stepping with " << tStepEnd << " time steps\n");

                wavefields->resetWavefields();

                ValueType DTinv = 1 / config.get<ValueType>("DT");

                start_t_shot = common::Walltime::get();
                for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {

                    *wavefieldsTemp = *wavefields;

                    solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);

                    // save wavefields in std::vector
                    *wavefieldrecord[tStep] = *wavefields;
                    //calculate temporal derivative of wavefield
                    *wavefieldrecord[tStep] -= *wavefieldsTemp;
                    *wavefieldrecord[tStep] *= DTinv;

                    if (config.get<bool>("useEnergyPreconditioning") == 1) {
                        energyPrecond.intSquaredWavefields(*wavefields, config.get<ValueType>("DT"));
                    }
                }

                receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration) + ".shot_" + std::to_string(shotNumber));

                HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << sources.getNumShots() << ": Calculate misfit and adjoint sources\n");

                /* Normalize observed and synthetic data */
                receivers.getSeismogramHandler().normalize();
                receiversTrue.getSeismogramHandler().normalize();

                /* Calculate misfit of one shot */
                misfitPerIt.setValue(shotNumber, dataMisfit->calc(receivers, receiversTrue));

                /* Calculate adjoint sources */
                dataMisfit->calcAdjointSources(adjointSources, receivers, receiversTrue);

                /* Calculate gradient */

                HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << sources.getNumShots() << ": Start Backward\n");
                gradientCalculation.run(*solver, *derivatives, receivers, sources, adjointSources, *model, *gradientPerShot, wavefieldrecord, config, shotNumber, workflow);

                /* Apply energy preconditioning per shot */
                if (config.get<bool>("useEnergyPreconditioning") == 1) {
                    energyPrecond.apply(*gradientPerShot, shotNumber);
                }

                if (config.get<bool>("useReceiversPerShot"))
                    ReceiverTaper.apply(*gradientPerShot);

                gradientPerShot->normalize();
                *gradient += *gradientPerShot;

                solver->resetCPML();

                end_t_shot = common::Walltime::get();
                HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << sources.getNumShots() << ": Finished in " << end_t_shot - start_t_shot << " sec.\n");

            } //end of loop over shots

            gradient->sumShotDomain(commInterShot);
            commInterShot->sumArray(misfitPerIt.getLocalValues());

            HOST_PRINT(commAll, "\n======== Finished loop over shots =========");
            HOST_PRINT(commAll, "\n===========================================\n");

            /* Apply receiver Taper (if ReceiverTaperRadius=0 gradient will be multplied by 1) */
            if (!config.get<bool>("useReceiversPerShot"))
                ReceiverTaper.apply(*gradient);

            gradientOptimization->apply(*gradient, workflow, *model);

            if (config.get<IndexType>("FreeSurface") == 2) {
                lama::DenseVector<ValueType> mask;
                mask = model->getVelocityP();
                mask.unaryOp(mask, common::UnaryOp::SIGN);
                *gradient *= mask;
            }

            if (config.get<bool>("useGradientTaper"))
                gradientTaper->apply(*gradient);

            /* Output of gradient */
            if (config.get<IndexType>("WriteGradient") && commInterShot->getRank() == 0)
                gradient->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1), config.get<IndexType>("PartitionedOut"), workflow);

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

            SLsearch.run(commAll, *solver, *derivatives, receivers, sources, receiversTrue, *model, dist, config, *gradient, steplengthInit, dataMisfit->getMisfitIt(workflow.iteration), workflow, freqFilter, sourceEst, *sourceSignalTaper);

            HOST_PRINT(commAll, "=========== Update Model ============\n\n");
            /* Apply model update */
            *gradient *= SLsearch.getSteplength();
            *model -= *gradient;

            if (config.get<bool>("useModelThresholds"))
                model->applyThresholds(config);

            if (commInterShot->getRank() == 0)
                model->write((config.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1)), config.get<IndexType>("PartitionedOut"));

            steplengthInit *= 0.98;

            end_t = common::Walltime::get();
            HOST_PRINT(commAll, "\nFinished iteration " << workflow.iteration + 1 << " in " << end_t - start_t << " sec.\n\n");

            /* -------------------------------------------------------------------- */
            /* One extra forward modelling to ensure complete and consistent output */
            /* -------------------------------------------------------------------- */
            if (workflow.iteration == maxiterations - 1) {

                HOST_PRINT(commAll, "================ Maximum number of iterations reached =================\n");
                HOST_PRINT(commAll, "=== Do one more forward modelling to calculate misfit and save seismograms ===\n\n");

                /* Update model for fd simulation (averaging, inverse Density ...) */
                model->prepareForModelling(config, ctx, dist, commShot);
                solver->prepareForModelling(*model, config.get<ValueType>("DT"));

                for (IndexType shotNumber = shotDist.lb(); shotNumber < shotDist.ub(); shotNumber++) {

                    /* Read field data (or pseudo-observed data, respectively) */
                    if (config.get<bool>("useReceiversPerShot")) {
                        receivers.init(config, ctx, dist, shotNumber);
                        receiversTrue.init(config, ctx, dist, shotNumber);
                    }
                    receiversTrue.getSeismogramHandler().read(config, fieldSeisName + ".shot_" + std::to_string(shotNumber), 1);

                    HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << sources.getNumShots() << ": Additional forward run with " << tStepEnd << " time steps\n");

                    wavefields->resetWavefields();

                    sources.init(config, ctx, dist, shotNumber);
                    if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0) {
                        sources.getSeismogramHandler().filter(freqFilter);
                        receiversTrue.getSeismogramHandler().filter(freqFilter);
                    }

                    if (config.get<bool>("useSourceSignalInversion") == 1) {
                        sourceEst.applyFilter(sources, shotNumber);
                        if (config.get<bool>("useSourceSignalTaper"))
                            sourceSignalTaper->apply(sources.getSeismogramHandler());
                    }

                    start_t_shot = common::Walltime::get();

                    for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                        solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);
                    }

                    receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber));

                    /* Normalize observed and synthetic data */
                    receivers.getSeismogramHandler().normalize();
                    receiversTrue.getSeismogramHandler().normalize();

                    /* Calculate misfit of one shot */
                    misfitPerIt.setValue(shotNumber, dataMisfit->calc(receivers, receiversTrue));

                    end_t_shot = common::Walltime::get();
                    HOST_PRINT(commShot, "Shot " << shotNumber + 1 << " of " << sources.getNumShots() << ": Finished additional forward run\n");

                } //end of loop over shots

                commInterShot->sumArray(misfitPerIt.getLocalValues());
                dataMisfit->addToStorage(misfitPerIt);

                SLsearch.appendToLogFile(commAll, workflow.workflowStage + 1, workflow.iteration + 1, logFilename, dataMisfit->getMisfitSum(workflow.iteration + 1));

                if (workflow.workflowStage != workflow.maxStage - 1) {
                    HOST_PRINT(commAll, "\nChange workflow stage\n");
                    workflow.changeStage(config, *dataMisfit, steplengthInit);
                }

            } // end extra forward modelling

        } // end of loop over iterations

    } // end of loop over workflow stages

    return 0;
}
