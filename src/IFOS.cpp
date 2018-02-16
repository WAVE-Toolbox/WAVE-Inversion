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
#include "Gradient/GradientFactory.hpp"
#include "Optimization/StepLengthSearch.hpp"

#include <Common/HostPrint.hpp>
#include <Partitioning/PartitioningCubes.hpp>

using namespace scai;
using namespace KITGPI;

int main(int argc, char *argv[])
{
    typedef double ValueType;
    double start_t, end_t; /* For timing */

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
    std::string fieldSeisName(config.get<std::string>("FieldSeisName"));
    std::string gradname(config.get<std::string>("GradientFilename"));
    std::string logFilename = config.get<std::string>("LogFilename");
    lama::Scalar steplength_init = config.get<ValueType>("SteplengthInit");

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
    KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefields->init(ctx, dist);
    
    /* --------------------------------------- */
    /* Wavefield record                        */
    /* --------------------------------------- */
    typedef typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
    std::vector<wavefieldPtr> wavefieldrecord;
    
    for (IndexType i = 0; i < tEnd; i++) {
        wavefieldPtr wavefieldsTemp(KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType));
        wavefieldsTemp->init(ctx, dist);

        wavefieldrecord.push_back(wavefieldsTemp);
    }
    
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
    gradientCalculation.allocate(config, dist, ctx);
    
    /* --------------------------------------- */
    /* Misfit                                  */
    /* --------------------------------------- */
    Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(Misfit::Factory<ValueType>::Create(misfitType));
    scai::lama::DenseVector<ValueType> misfitPerIt(sources.getNumShots(), 0, ctx);
    
    /* --------------------------------------- */
    /* Step length search                      */
    /* --------------------------------------- */
    StepLengthSearch<ValueType> SLsearch;
    SLsearch.initLogFile(comm, logFilename);

    /* --------------------------------------- */
    /* Gradient taper                          */
    /* --------------------------------------- */
    SourceReceiverTaper<ValueType> ReceiverTaper;
    ReceiverTaper.init(dist,ctx,receivers,config,config.get<IndexType>("SourceTaperRadius"));

    /* --------------------------------------- */
    /* True data                               */
    /* --------------------------------------- */
    KITGPI::Acquisition::Receivers<ValueType> receiversTrue(config, ctx, dist);
    
    
    /* --------------------------------------- */
    /*        Loop over iterations             */
    /* --------------------------------------- */

    IndexType maxiterations = config.get<IndexType>("MaxIterations");
    if (config.get<bool>("runForward"))
        maxiterations = 1;
    for (IndexType iteration = 0; iteration < maxiterations; iteration++) {

        // update model for fd simulation (averaging, inverse Density ...)
        model->prepareForModelling(config, ctx, dist, comm);
        
        /* --------------------------------------- */
        /*        Loop over shots                  */
        /* --------------------------------------- */
        
        gradient->reset();    // reset gradient because gradient is a sum of all gradientsPerShot gradients+=gradientPerShot

        for (IndexType shotNumber = 0; shotNumber < sources.getNumShots(); shotNumber++) {
                   
            /* Read field data (or pseudo-observed data, respectively) */
            receiversTrue.getSeismogramHandler().readFromFileRaw(fieldSeisName + ".It0" + ".shot" + std::to_string(shotNumber) + ".mtx", 1);
            
            /* --------------------------------------- */
            /*        Forward modelling                */
            /* --------------------------------------- */
            
            HOST_PRINT(comm, "\n================Start Forward====================\n");
            HOST_PRINT(comm, "Start time stepping for shot " << shotNumber + 1 << " of " << sources.getNumShots() << "\n"
                                                                << "Total Number of time steps: " << tEnd << "\n");
            wavefields->reset();
            sources.init(config, ctx, dist, shotNumber);

            start_t = scai::common::Walltime::get();
            for (t = 0; t < tEnd; t++) {

                solver->run(receivers, sources, *model, *wavefields, *derivatives, t, t + 1, config.get<ValueType>("DT"));

                // save wavefields in Dense Matrix
                wavefieldrecord[t]->assign(*wavefields);
            }

            receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename") + ".It" + std::to_string(iteration) + ".shot" + std::to_string(shotNumber));

            end_t = scai::common::Walltime::get();
            HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");
            
            /* Calculate misfit of one shot */
            misfitPerIt.setValue(shotNumber, dataMisfit->calc(receivers, receiversTrue));          
            
            /* Calculate adjoint sources */
 
            /* Calculate gradient */
            gradientCalculation.calc(*solver, *derivatives, receivers, sources, receiversTrue, *model, *gradientPerShot, wavefieldrecord, config, iteration, shotNumber);
            *gradient += *gradientPerShot; 
        
        } //end of loop over shots
        
        dataMisfit->add(misfitPerIt);
        
        HOST_PRINT(comm, "Misfit after iteration " << iteration << ": " << dataMisfit->getMisfitSum(iteration) << "\n\n");

        if ((iteration > 0) && (dataMisfit->getMisfitSum(iteration) >= dataMisfit->getMisfitSum(iteration - 1))) {
            HOST_PRINT(comm, "Misfit is getting higher after iteration " << iteration << ", last_misfit: " << dataMisfit->getMisfitSum(iteration - 1) << "\n\n");
            break;
        }
        
        /* Apply receiver Taper (if ReceiverTaperRadius=0 gradient will be multplied by 1) */
        ReceiverTaper.apply(*gradient);
        
        /* Output of gradient */
        if(config.get<IndexType>("WriteGradient"))
            gradient->getVelocityP().writeToFile(gradname + "_vp" + ".It" + std::to_string(iteration) + ".mtx");
       
        gradient->scale(*model);
        
        SLsearch.calc(*solver, *derivatives, receivers, sources, receiversTrue, *model, dist, config, *gradient, steplength_init, dataMisfit->getMisfitSum(iteration));
        
        *gradient *= SLsearch.getSteplength();
        *model -= *gradient;
       
        model->write((config.get<std::string>("ModelFilename") + ".It" + std::to_string(iteration)), config.get<IndexType>("PartitionedOut"));

        steplength_init*=0.98; // 0.95 with steplengthMax = 0.1 yields misfit of ~97 
        
        SLsearch.appendToLogFile(comm, iteration, logFilename);
 
    } //end of loop over iterations
    return 0;
}
