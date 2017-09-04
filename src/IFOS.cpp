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

    std::string convname("gradients/waveconv");
    std::string gradname("gradients/grad");

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

    //     /* --------------------------------------- */
    //     /* Wavefields                              */
    //     /* --------------------------------------- */
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(Wavefields::Factory<ValueType>::Create(dimension, equationType));
    wavefields->init(ctx, dist);
    //
    //     /* --------------------------------------- */
    //     /* Acquisition geometry                    */
    //     /* --------------------------------------- */
    Acquisition::Receivers<ValueType> receivers(config, ctx, dist);
    Acquisition::Sources<ValueType> sources(config, ctx, dist);
    //  sources.init;

    //     /* --------------------------------------- */
    //     /* Modelparameter                          */
    //     /* --------------------------------------- */
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    model->init(config, ctx, dist);
    lama::Vector &VelocityP = model->getVelocityP();
  //  lama::Vector &Density = model->getDensity();

    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */
    IndexType getNT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));
    solver->prepareBoundaryConditions(config, *derivatives, dist, ctx);

    /* --------------------------------------- */
    /* Wavefields and Convulution fiels        */
    /* --------------------------------------- */

    lama::DenseVector<ValueType> waveconv_v;
    lama::DenseVector<ValueType> waveconv_p;
    lama::DenseVector<ValueType> tmp;
    lama::DenseVector<ValueType> waveconv_v_sum;
    lama::DenseVector<ValueType> waveconv_p_sum;

    waveconv_v.allocate(dist);
    waveconv_v.setContextPtr(ctx);
    waveconv_p.allocate(dist);
    waveconv_p.setContextPtr(ctx);
    tmp.allocate(dist);
    tmp.setContextPtr(ctx);
    waveconv_v_sum.allocate(dist);
    waveconv_v_sum.setContextPtr(ctx);
    waveconv_p_sum.allocate(dist);
    waveconv_p_sum.setContextPtr(ctx);

    lama::DenseMatrix<ValueType> wavefieldrecordvx;
    lama::DenseMatrix<ValueType> wavefieldrecordvy;
    lama::DenseMatrix<ValueType> wavefieldrecordp;

    dmemo::DistributionPtr no_dist_NT(new scai::dmemo::NoDistribution(getNT));

    wavefieldrecordvx.allocate(dist, no_dist_NT);
    wavefieldrecordvx.setContextPtr(ctx);
    wavefieldrecordvy.allocate(dist, no_dist_NT);
    wavefieldrecordvy.setContextPtr(ctx);
    wavefieldrecordp.allocate(dist, no_dist_NT);
    wavefieldrecordp.setContextPtr(ctx);

    lama::Scalar Misfit;
    lama::Scalar MisfitSum;
    lama::Scalar MisfitSumTemp;
    
    ValueType steplength=80;
    
    /* --------------------------------------- */
    /*               Gradients                 */
    /* --------------------------------------- */

    lama::DenseVector<ValueType> grad_bulk;
    lama::DenseVector<ValueType> grad_rho;
    lama::DenseVector<ValueType> grad_vp;

    grad_vp.allocate(dist);
    grad_vp.setContextPtr(ctx);
    grad_bulk.allocate(dist);
    grad_bulk.setContextPtr(ctx);
    grad_rho.allocate(dist);
    grad_rho.setContextPtr(ctx);

    /* --------------------------------------- */
    /* Adjoint sources:                         */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> adjoint(config, ctx, dist);
    Acquisition::Seismogram<ValueType> truedata(adjoint.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P));
    Acquisition::Seismogram<ValueType> synthetic(receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P));

    // =================================== Forward Modelling ==============================================

    /* Start and end counter for time stepping */
    IndexType t = 0;
    IndexType tEnd = getNT;
    IndexType maxiterations=10;
if (config.get<bool>("runForward"))
    maxiterations=1;
    for (IndexType iteration = 0; iteration < maxiterations; iteration++) {
             MisfitSumTemp=0;
	     model->prepareForModelling(config, ctx, dist, comm);
	     waveconv_p_sum.assign(0);
	     
	     
        for (IndexType shotNumber = 0; shotNumber < sources.getNumShots(); shotNumber++) {
            /* Update Source */
            if (!config.get<bool>("runSimultaneousShots")) {
                sources.init(config, ctx, dist, shotNumber);
            }

            wavefields->reset();
            waveconv_p = 0;
            waveconv_v = 0;

            HOST_PRINT(comm, "\n================Start Forward====================\n");
            HOST_PRINT(comm, "Start time stepping for shot " << shotNumber + 1 << " of " << sources.getNumShots() << "\n"
                                                             << "Total Number of time steps: " << getNT << "\n");
            start_t = common::Walltime::get();
            for (t = 0; t < tEnd; t++) {

                solver->run(receivers, sources, *model, *wavefields, *derivatives, t, t + 1, config.get<ValueType>("DT"));

                // save wavefields in Dense Matrix
                wavefieldrecordvx.setColumn(wavefields->getVX(), t, scai::common::binary::BinaryOp::COPY);
                wavefieldrecordvy.setColumn(wavefields->getVY(), t, scai::common::binary::BinaryOp::COPY);
                wavefieldrecordp.setColumn(wavefields->getP(), t, scai::common::binary::BinaryOp::COPY);
            }
            end_t = common::Walltime::get();
            HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");

            if (!config.get<bool>("runSimultaneousShots")) {
                receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename") +".It" + std::to_string(iteration) +".shot" + std::to_string(shotNumber));
            } else {
                receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename"));
                //wavefieldrecord.writeToFile("wavefields/test.mtx");
            }
            // =========================================================================================================

            if (!config.get<bool>("runForward")) {

                /* --------------------------------------- */
                /* Adjoint sources:                         */
                /* --------------------------------------- */
                std::string FiledSeisName("seismograms/rectangle.true");
                truedata.readFromFileRaw(FiledSeisName +".It0" + ".shot" + std::to_string(shotNumber) + ".p.mtx", adjoint.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData().getRowDistributionPtr(), NULL);
                synthetic = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P);


                synthetic -= truedata;
		Misfit=0.5*synthetic.getData().l2Norm();
		
                adjoint.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P) = synthetic;
                //  adjoint.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).writeToFileRaw("seismograms/adjoint.mtx");

                //======================================Backward Modelling==================================

                wavefields->reset();

                HOST_PRINT(comm, "\n================Start Backward====================\n");
                HOST_PRINT(comm, "Start time stepping\n"
                                     << "Total Number of time steps: " << getNT << "\n");
                start_t = common::Walltime::get();

                for (t = tEnd - 1; t >= 0; t--) {

                    solver->run(receivers, adjoint, *model, *wavefields, *derivatives, t, t + 1, config.get<ValueType>("DT"));

                    /* --------------------------------------- */
                    /* Convolution:                         */
                    /* --------------------------------------- */

                    wavefieldrecordp.getColumn(tmp, t);
                    tmp *= wavefields->getP();
                    waveconv_p += tmp;

//                     wavefieldrecordvx.getColumn(tmp, t);
//                     tmp *= wavefields->getVX();
//                     waveconv_v += tmp;
//                     wavefieldrecordvy.getColumn(tmp, t);
//                     tmp *= wavefields->getVY();
//                     waveconv_v += tmp;
                }
                
                MisfitSumTemp+=Misfit;
                waveconv_p_sum += waveconv_p;
        //        waveconv_v_sum += waveconv_v;

                end_t = common::Walltime::get();
                HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");

                //   receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).writeToFileRaw("seismograms/rec_adjoint.mtx");
                //
               //     waveconv_p.writeToFile(convname + ".shot_" + std::to_string(shotNumber)+".mtx");
                //=========================================================================================================
            }
        } //end of loop over shots
        
              
        if (!config.get<bool>("runForward")) {
		HOST_PRINT(comm,"Misfit " << MisfitSumTemp << "iteration " << iteration << "\n\n" );
	  
		if ((iteration>0) && (MisfitSumTemp>MisfitSum))
		{
			HOST_PRINT(comm,"Misfit is getting higher after iteration: " << iteration << "last_misfit: " << MisfitSum <<"\n\n" );
		break;	
		}
		MisfitSum=MisfitSumTemp;
            //    Output jacobi
            waveconv_p_sum.writeToFile(convname + "_p" +".It"+std::to_string(iteration)+ ".mtx");
//             waveconv_v_sum.writeToFile(convname + "_v" + ".mtx");

            //calculate gradient vp
            grad_bulk = model->getPWaveModulus();
            grad_bulk *= grad_bulk;
            grad_bulk.invert();
            grad_bulk *= waveconv_p_sum;
            grad_bulk *= -config.get<ValueType>("DT");

            grad_vp = 2 * grad_bulk;
            grad_vp *= model->getDensity();
            grad_vp *= model->getVelocityP();

//             grad_rho = model->getVelocityP();
//             grad_rho *= model->getVelocityP();
//             grad_rho *= grad_bulk;
//             grad_rho -= config.get<ValueType>("DT") * waveconv_v_sum;

            grad_vp.writeToFile(gradname + "_vp" + ".It"+std::to_string(iteration)+ ".mtx");
    //        grad_rho.writeToFile(gradname + "_rho" + ".It"+std::to_string(iteration)+ ".mtx");

	    steplength*=0.8;
            grad_vp *= 1 / grad_vp.max() * (steplength);
      //      grad_rho *= 1 / grad_rho.max() * 50;

            //  temp-=grad_vp;
            VelocityP -= grad_vp;
          //  Density -= grad_rho;

            //  model->getDensity()=grad_rho;

            model->write((config.get<std::string>("ModelFilename") + ".It"+std::to_string(iteration)), config.get<IndexType>("PartitionedOut"));
	    
        }
    } //end of loop over iterations
    return 0;
}
