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

#include "Optimization/GradientCalculation.hpp"

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

    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */
    /* inter node communicator */
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(comm->getNodeRank());
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    /* inter node distribution */
    // block distribution: i-st processor gets lines [i * N/num_processes] to [(i+1) * N/num_processes - 1] of the matrix
    IndexType getN = config.get<IndexType>("NZ") * config.get<IndexType>("NX") * config.get<IndexType>("NY");
    dmemo::DistributionPtr dist(new dmemo::BlockDistribution(getN, comm));

    if (config.get<IndexType>("UseCubePartitioning")) {
        Partitioning::PartitioningCubes<ValueType> partitioning(config, comm);
        dist = partitioning.getDist();
    }

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
    /* Wavefields                              */
    /* --------------------------------------- */
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(Wavefields::Factory<ValueType>::Create(dimension, equationType));
    wavefields->init(ctx, dist);

    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> receivers(config, ctx, dist);
    Acquisition::Sources<ValueType> sources(config, ctx, dist);
    //  sources.init;

    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    model->init(config, ctx, dist);

    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */
    IndexType getNT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));
    solver->prepareBoundaryConditions(config, *derivatives, dist, ctx);

    dmemo::DistributionPtr no_dist_NT(new scai::dmemo::NoDistribution(getNT));

//     lama::GridVector<ValueType> wavefieldStorage;
        
    /* --------------------------------------- */
    /* Model update                            */
    /* --------------------------------------- */

    GradientCalculation<ValueType> gradient;                  
    gradient.allocate(dist, no_dist_NT, comm, ctx);
    
    lama::DenseVector<ValueType> modelupdate(dist,ctx);
    
    ValueType steplength=80;
    
    /* --------------------------------------- */
    /*        Loop over iterations             */
    /* --------------------------------------- */
    
    IndexType maxiterations=10;
    if (config.get<bool>("runForward"))
        maxiterations=1;
    for (IndexType iteration = 0; iteration < maxiterations; iteration++) {
        
         model->prepareForModelling(config, ctx, dist, comm);  
         gradient.calc(*solver, *derivatives, receivers, sources, *model, *wavefields, config, iteration);
         
// 		HOST_PRINT(comm,"Misfit " << MisfitSumTemp << "iteration " << iteration << "\n\n" );
// 	  
// 		if ((iteration>0) && (MisfitSumTemp>MisfitSum))
// 		{
// 			HOST_PRINT(comm,"Misfit is getting higher after iteration: " << iteration << "last_misfit: " << MisfitSum <<"\n\n" );
// 		break;	
// 		}
          
	    steplength*=0.8;
        gradient.grad_vp *= 1 / gradient.grad_vp.max() * (steplength);
//        grad_rho *= 1 / grad_rho.max() * 50;

        modelupdate=model->getVelocityP()-gradient.grad_vp;
        //  temp-=grad_vp;
        model->setVelocityP(modelupdate);
        //  Density -= grad_rho;

        //  model->getDensity()=grad_rho;

        model->write((config.get<std::string>("ModelFilename") + ".It"+std::to_string(iteration)), config.get<IndexType>("PartitionedOut"));
	    
    } //end of loop over iterations
    return 0;
}