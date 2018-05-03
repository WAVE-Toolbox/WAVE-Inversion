#include <scai/lama.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include <Configuration/Configuration.hpp>
//#include "/home/tmetz/projects/FDSimulation_LAMA//src/Configuration/Configuration.hpp"
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>

#include <ForwardSolver/ForwardSolver.hpp>

#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolver/ForwardSolverFactory.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>

#include <Wavefields/WavefieldsFactory.hpp>

#include "Gradient/GradientFactory.hpp"
#include "Optimization/GradientCalculation.hpp"
#include "Optimization/Misfit/Misfit.hpp"
#include "Optimization/Misfit/MisfitFactory.hpp"
#include "Optimization/StepLengthSearch.hpp"
#include "Optimization/Preconditioning/EnergyPreconditioning.hpp"
#include "Workflow/Workflow.hpp"

#include <Common/HostPrint.hpp>
#include <Partitioning/PartitioningCubes.hpp>

using namespace scai;
using namespace KITGPI;

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
    IndexType t = 0;
    IndexType tEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
    
    std::string misfitType = config.get<std::string>("misfitType");
    std::string fieldSeisName(config.get<std::string>("fieldSeisName"));
    std::string gradname(config.get<std::string>("gradientFilename"));
    std::string logFilename = config.get<std::string>("logFilename");
    ValueType steplengthInit = config.get<ValueType>("steplengthInit");
    IndexType maxiterations = config.get<IndexType>("maxIterations");

    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */
    /* inter node communicator */
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(comm->getNodeRank());
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    // inter node distribution
    // define the grid topology by sizes NX, NY, and NZ from configuration
    // Attention: LAMA uses row-major indexing while SOFI-3D uses column-major, so switch dimensions, x-dimension has stride 1

    common::Grid3D grid(config.get<IndexType>("NZ"), config.get<IndexType>("NY"), config.get<IndexType>("NX"));
    common::Grid3D procGrid(config.get<IndexType>("ProcNZ"), config.get<IndexType>("ProcNY"), config.get<IndexType>("ProcNX"));
    // distribute the grid onto available processors, topology can be set by environment variable
    dmemo::DistributionPtr dist(new dmemo::GridDistribution(grid, comm, procGrid));

    HOST_PRINT(comm, "\nIFOS" << dimension << " " << equationType << " - LAMA Version\n\n");
    if (comm->getRank() == MASTERGPI) {
        config.print();
    }

    /* --------------------------------------- */
    /* Calculate derivative matrices           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    derivatives->init(dist, ctx, config, comm);
    end_t = common::Walltime::get();
    HOST_PRINT(comm, "Finished initializing matrices in " << end_t - start_t << " sec.\n\n");

    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> receivers(config, ctx, dist);
    Acquisition::Sources<ValueType> sources(config, ctx, dist);

    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    // load starting model
    model->init(config, ctx, dist);
    
    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */
    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));
    solver->prepareBoundaryConditions(config, *derivatives, dist, ctx);
    
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
    
    for (IndexType i = 0; i < tEnd; i++) {
        wavefieldPtr wavefieldsTemp(Wavefields::Factory<ValueType>::Create(dimension, equationType));
        wavefieldsTemp->init(ctx, dist);

        wavefieldrecord.push_back(wavefieldsTemp);
    }
    
    /* --------------------------------------- */
    /* True data                               */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> receiversTrue(config, ctx, dist);
    
    /* --------------------------------------- */
    /* Misfit                                  */
    /* --------------------------------------- */
    Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(Misfit::Factory<ValueType>::Create(misfitType));
    lama::DenseVector<ValueType> misfitPerIt(sources.getNumShots(), 0, ctx);
    
    /* --------------------------------------- */
    /* Adjoint sources                         */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> adjointSources(config, ctx, dist);
    
    /* --------------------------------------- */
    /* Workflow                                */
    /* --------------------------------------- */
    Workflow::Workflow<ValueType> workflow(config);
    
    /* --------------------------------------- */
    /* Step length search                      */
    /* --------------------------------------- */
    StepLengthSearch<ValueType> SLsearch;
    SLsearch.initLogFile(comm, logFilename, misfitType);
    
    /* --------------------------------------- */
    /* Gradients                               */
    /* --------------------------------------- */
    Gradient::Gradient<ValueType>::GradientPtr gradient(Gradient::Factory<ValueType>::Create(equationType));
    gradient->init(ctx, dist);
    
    Gradient::Gradient<ValueType>::GradientPtr gradientPerShot(Gradient::Factory<ValueType>::Create(equationType));
    gradientPerShot->init(ctx, dist);
    
    /* --------------------------------------- */
    /* Gradient calculation                    */
    /* --------------------------------------- */
    GradientCalculation<ValueType> gradientCalculation;
    
    /* --------------------------------------- */
    /* Gradient taper                          */
    /* --------------------------------------- */
    Preconditioning::SourceReceiverTaper<ValueType> ReceiverTaper;
    ReceiverTaper.init(dist,ctx,receivers,config,config.get<IndexType>("SourceTaperRadius"));
    
    /* --------------------------------------- */
    /* Gradient preconditioning                */
    /* --------------------------------------- */
    Preconditioning::EnergyPreconditioning<ValueType> energyPrecond;
    if (config.get<bool>("useEnergyPreconditioning") == 1){
        energyPrecond.init(dist, config);}
    
    
    /* --------------------------------------- */
    /*       Loop over workflow stages         */
    /* --------------------------------------- */
    
    for (IndexType workflowStage = 0; workflowStage < workflow.maxStage; workflowStage++) {
        
        HOST_PRINT(comm, "\n=================================================");
        HOST_PRINT(comm, "\n============ Workflow stage " << workflowStage+1 << " of " << workflow.maxStage << " ==============");    
        HOST_PRINT(comm, "\n=================================================\n\n");
        
        HOST_PRINT(comm, "Set parameters: \n");
        HOST_PRINT(comm, "invertForVp = " << workflow.invertForVp << "\n");
        HOST_PRINT(comm, "invertForVs = " << workflow.invertForVs << "\n");
        HOST_PRINT(comm, "invertForDensity = " << workflow.invertForDensity << "\n");
    
        gradientCalculation.allocate(config, dist, ctx, workflow); 
        
        /* --------------------------------------- */
        /*        Loop over iterations             */
        /* --------------------------------------- */

        for (IndexType iteration = 0; iteration < maxiterations; iteration++) {
            
            HOST_PRINT(comm, "\n=================================================");
            HOST_PRINT(comm, "\n================ Iteration " << iteration+1 << " ====================");    
            HOST_PRINT(comm, "\n=================================================\n\n");
            start_t = common::Walltime::get();
            
            /* Update model for fd simulation (averaging, inverse Density ...) */
            model->prepareForModelling(config, ctx, dist, comm);
            
            /* --------------------------------------- */
            /*        Loop over shots                  */
            /* --------------------------------------- */
            
            gradient->resetGradient();    // reset gradient because gradient is a sum of all gradientsPerShot gradients+=gradientPerShot

            for (IndexType shotNumber = 0; shotNumber < sources.getNumShots(); shotNumber++) {
                
                gradientPerShot->resetGradient(); 
                    
                /* Read field data (or pseudo-observed data, respectively) */
                receiversTrue.getSeismogramHandler().readFromFileRaw(fieldSeisName + ".shot_" + std::to_string(shotNumber) + ".mtx", 1);
                
                /* Reset approximated Hessian per shot */
                if (config.get<bool>("useEnergyPreconditioning") == 1){
                    energyPrecond.resetApproxHessian();}
                
                HOST_PRINT(comm, "\n=============== Shot " << shotNumber + 1 << " of " << sources.getNumShots() << " ===================\n");
            
                /* --------------------------------------- */
                /*        Forward modelling                */
                /* --------------------------------------- */
                

                HOST_PRINT(comm, "\n================Start Forward====================\n");
                HOST_PRINT(comm, "Start time stepping for shot " << shotNumber + 1 << " of " << sources.getNumShots() << "\n"
                                                                    << "Total Number of time steps: " << tEnd << "\n");
                                                                    
                wavefields->resetWavefields();

                sources.init(config, ctx, dist, shotNumber);
            
                ValueType DTinv=1/config.get<ValueType>("DT");
            
                start_t_shot = common::Walltime::get();
                for (t = 0; t < tEnd; t++) {
                    
                    *wavefieldsTemp=*wavefields;
                    
                    solver->run(receivers, sources, *model, *wavefields, *derivatives, t, t + 1, config.get<ValueType>("DT"));

                    // save wavefields in std::vector
                    *wavefieldrecord[t]=*wavefields;
                    //calculate temporal derivative of wavefield
                    *wavefieldrecord[t]-=*wavefieldsTemp;
                    *wavefieldrecord[t]*=DTinv;
                    
                    if (config.get<bool>("useEnergyPreconditioning") == 1){
                        energyPrecond.intSquaredWavefields(*wavefields, config.get<ValueType>("DT"));}
                }

                receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflowStage+1) + ".It_" + std::to_string(iteration) + ".shot_" + std::to_string(shotNumber));
                
                HOST_PRINT(comm, "\nCalculate misfit and adjoint sources\n");
                
                /* Calculate misfit of one shot */
                misfitPerIt.setValue(shotNumber, dataMisfit->calc(receivers, receiversTrue));
                
                /* Calculate adjoint sources */
                dataMisfit->calcAdjointSources(adjointSources, receivers, receiversTrue);
                
                /* Calculate gradient */
                gradientCalculation.run(*solver, *derivatives, receivers, sources, adjointSources, *model, *gradientPerShot, wavefieldrecord, config, workflowStage, iteration, shotNumber, workflow);
                
                /* Apply energy preconditioning per shot */
                if (config.get<bool>("useEnergyPreconditioning") == 1){
                    energyPrecond.apply(*gradientPerShot, shotNumber);}
                
                *gradient += *gradientPerShot; 
                
                end_t_shot = common::Walltime::get();
                HOST_PRINT(comm, "\nFinished shot in " << end_t_shot - start_t_shot << " sec.\n\n");
            
            } //end of loop over shots
            
            HOST_PRINT(comm, "\n======== Finished loop over shots =========");
            HOST_PRINT(comm, "\n===========================================\n");

            /* Output of gradient */
            if(config.get<IndexType>("WriteGradient"))
                gradient->write(gradname + ".stage_" + std::to_string(workflowStage+1) + ".It_" + std::to_string(iteration + 1), config.get<IndexType>("PartitionedOut"));

            dataMisfit->addToStorage(misfitPerIt);
            
            SLsearch.appendToLogFile(comm, workflowStage+1, iteration, logFilename, dataMisfit->getMisfitSum(iteration));
            
            /* Check abort criteria */
            HOST_PRINT(comm, "\nMisfit after stage " << workflowStage << ", iteration " << iteration << ": " << dataMisfit->getMisfitSum(iteration) << "\n");

            if ((iteration > 0) && (dataMisfit->getMisfitSum(iteration) >= dataMisfit->getMisfitSum(iteration - 1))) {
                HOST_PRINT(comm, "\nMisfit is getting higher after iteration " << iteration << ", last_misfit: " << dataMisfit->getMisfitSum(iteration - 1) << "\n\n");
                workflow.changeStage(config, *dataMisfit, steplengthInit);
                break;
            }
            
            /* Apply receiver Taper (if ReceiverTaperRadius=0 gradient will be multplied by 1) */
            ReceiverTaper.apply(*gradient);
            
            gradient->scale(*model, workflow);
            
            HOST_PRINT(comm,"\n===========================================" );
            HOST_PRINT(comm,"\n======== Start step length search =========\n" );
        
            SLsearch.run(*solver, *derivatives, receivers, sources, receiversTrue, *model, dist, config, *gradient, steplengthInit, dataMisfit->getMisfitIt(iteration));
            
        
            HOST_PRINT(comm, "=========== Update Model ============\n\n");
            /* Apply model update */
            *gradient *= SLsearch.getSteplength();
            *model -= *gradient;
        
            model->write((config.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflowStage+1) + ".It_" + std::to_string(iteration+1)), config.get<IndexType>("PartitionedOut"));

            steplengthInit*=0.98; 
        
            end_t = common::Walltime::get();
            HOST_PRINT(comm, "\nFinished iteration " << iteration + 1 << " in " << end_t - start_t << " sec.\n\n");
            
            
            /* -------------------------------------------------------------------- */
            /* One extra forward modelling to ensure complete and consistent output */
            /* -------------------------------------------------------------------- */
            if (iteration == maxiterations-1){
                
                HOST_PRINT(comm, "================ Maximum number of iterations reached =================\n");
                HOST_PRINT(comm, "=== Do one more forward modelling to calculate misfit and save seismograms ===\n\n");
                
                /* Update model for fd simulation (averaging, inverse Density ...) */
                model->prepareForModelling(config, ctx, dist, comm);
                
                for (IndexType shotNumber = 0; shotNumber < sources.getNumShots(); shotNumber++) {
                    
                /* Read field data (or pseudo-observed data, respectively) */
                receiversTrue.getSeismogramHandler().readFromFileRaw(fieldSeisName + ".shot_" + std::to_string(shotNumber) + ".mtx", 1);  

                HOST_PRINT(comm, "\n================Start Forward====================\n");
                HOST_PRINT(comm, "Start time stepping for shot " << shotNumber + 1 << " of " << sources.getNumShots() << "\n"
                                                                    << "Total Number of time steps: " << tEnd << "\n");
                                                                    
                wavefields->resetWavefields();

                sources.init(config, ctx, dist, shotNumber);
            
                start_t_shot = common::Walltime::get();
 
                solver->run(receivers, sources, *model, *wavefields, *derivatives, 0, tEnd, config.get<ValueType>("DT"));

                receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflowStage+1) + ".It_" + std::to_string(iteration+1) + ".shot_" + std::to_string(shotNumber));
                
                /* Calculate misfit of one shot */
                misfitPerIt.setValue(shotNumber, dataMisfit->calc(receivers, receiversTrue));          
                
                end_t_shot = common::Walltime::get();
                HOST_PRINT(comm, "\nFinished shot in " << end_t_shot - start_t_shot << " sec.\n\n");
                
                } //end of loop over shots
                
                dataMisfit->addToStorage(misfitPerIt);
                
                SLsearch.appendToLogFile(comm, workflowStage+1, iteration+1, logFilename, dataMisfit->getMisfitSum(iteration+1));
                
                if(workflowStage != workflow.maxStage-1){
                    workflow.changeStage(config, *dataMisfit, steplengthInit);}
                
            } // end extra forward modelling
            
        } // end of loop over iterations
    
    } // end of loop over workflow stages
    
    return 0;
}
