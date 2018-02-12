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
    std::string misfitType = config.get<std::string>("misfitType");

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
    /* Calculate derivative matrizes           */
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
    /* Objects for inversion                   */
    /* --------------------------------------- */

    Gradient::Gradient<ValueType>::GradientPtr gradient(Gradient::Factory<ValueType>::Create(equationType));
    gradient->init(ctx, dist);

    GradientCalculation<ValueType> gradientCalculation;
    Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(Misfit::Factory<ValueType>::Create(misfitType));
    StepLengthSearch<ValueType> SLsearch;
    std::string logFilename = config.get<std::string>("LogFilename");
    SLsearch.initLogFile(comm, logFilename);

    gradientCalculation.allocate(config, dist, ctx);
    SourceReceiverTaper<ValueType> ReceiverTaper;
    ReceiverTaper.init(dist,ctx,receivers,config,config.get<IndexType>("SourceTaperRadius"));


    lama::Scalar steplength_init = config.get<ValueType>("SteplengthInit");
    /* --------------------------------------- */
    /*        Loop over iterations             */
    /* --------------------------------------- */
    std::string gradname(config.get<std::string>("GradientFilename"));

    IndexType maxiterations = config.get<IndexType>("MaxIterations");
    if (config.get<bool>("runForward"))
        maxiterations = 1;
    for (IndexType iteration = 0; iteration < maxiterations; iteration++) {

        // update model for fd simulation (averaging, inverse Density ...
        model->prepareForModelling(config, ctx, dist, comm);

        gradientCalculation.calc(*solver, *derivatives, receivers, sources, *model, *gradient, config, iteration, *dataMisfit);

        HOST_PRINT(comm, "Misfit after iteration " << iteration << ": " << dataMisfit->getMisfitSum(iteration) << "\n\n");

        //         abortcriterion.check(dataMisfit);
        if ((iteration > 0) && (dataMisfit->getMisfitSum(iteration) >= dataMisfit->getMisfitSum(iteration - 1))) {
            HOST_PRINT(comm, "Misfit is getting higher after iteration " << iteration << ", last_misfit: " << dataMisfit->getMisfitSum(iteration - 1) << "\n\n");
            break;
        }

        //apply receiver Taper (if ReceiverTaperRadius=0 gradient will be multplied by 1)
        ReceiverTaper.apply(*gradient);
	gradient->getVelocityP().writeToFile(gradname + "_vp" + ".It" + std::to_string(iteration) + ".mtx");
	 gradient->scale(*model);
        
        SLsearch.calc(*solver, *derivatives, receivers, sources, *model, dist, config, *gradient, steplength_init, dataMisfit->getMisfitSum(iteration));
        
        *gradient *=SLsearch.getSteplength();
        *model -= *gradient;
       
	
        model->write((config.get<std::string>("ModelFilename") + ".It" + std::to_string(iteration)), config.get<IndexType>("PartitionedOut"));

        steplength_init*=0.98; // 0.95 with steplengthMax = 0.1 yields misfit of ~97 
        
        SLsearch.appendToLogFile(comm, iteration, logFilename);
 
    } //end of loop over iterations
    return 0;
}
