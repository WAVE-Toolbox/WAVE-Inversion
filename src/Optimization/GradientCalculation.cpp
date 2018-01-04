#include <iostream>

#include "GradientCalculation.hpp"

template <typename ValueType>
void GradientCalculation<ValueType>::allocate(KITGPI::Configuration::Configuration config,scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr no_dist_NT, scai::dmemo::CommunicatorPtr comm, scai::hmemo::ContextPtr ctx)
{
        
	std::string dimension = config.get<std::string>("dimension");
        std::string equationType = config.get<std::string>("equationType");
	IndexType getNT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
	
	wavefields=KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
	wavefields->init(ctx, dist);
	
    /* ------------------------------------------------------- */
    /* Allocate wavefield record                               */
    /* ------------------------------------------------------- */
	
        for (IndexType i=0;i<getNT;i++){
        wavefieldPtr wavefieldsTemp(KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType));
        wavefieldsTemp->init(ctx, dist);
    
        wavefieldrecord.push_back(wavefieldsTemp);
    }
	
    /* ------------------------------------------- */
    /* Allocate convolution fields and set context */
    /* ------------------------------------------- */

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
    
    /* ----------------------------------------- */
    /* Set distribution and context of gradients */
    /* ----------------------------------------- */
    
    grad_bulk.allocate(dist);
    grad_bulk.setContextPtr(ctx);
    grad_rho.allocate(dist);
    grad_rho.setContextPtr(ctx);
    grad_vp.allocate(dist);
    grad_vp.setContextPtr(ctx);
    

}

template <typename ValueType>
void GradientCalculation<ValueType>::calc(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Configuration::Configuration config, IndexType iteration, Misfit<ValueType> &dataMisfit)
{
    
    double start_t, end_t; /* For timing */
    IndexType t = 0;
    IndexType getNT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);  
    
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
       
    scai::dmemo::DistributionPtr dist = wavefields->getVX().getDistributionPtr();
    scai::dmemo::DistributionPtr no_dist_NT(new scai::dmemo::NoDistribution(getNT));
    scai::dmemo::CommunicatorPtr comm = scai::dmemo::Communicator::getCommunicatorPtr();   // default communicator, set by environment variable SCAI_COMMUNICATOR
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                   // default context, set by environment variable SCAI_CONTEXT   

    std::string convname("gradients/waveconv");
    std::string gradname("gradients/grad");
    
    /* Should misfit be part of #GradientCalculation in a later stage? */
    scai::lama::DenseVector<ValueType> misfitTemp(sources.getNumShots(), 0, ctx);

    /* ------------------------------------------------------ */
    /* Allocate adjoint sources, true data and synthetic data */
    /* ------------------------------------------------------ */
    
    KITGPI::Acquisition::Receivers<ValueType> adjoint(config, ctx, dist);
    KITGPI::Acquisition::Seismogram<ValueType> truedata(adjoint.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P));
    KITGPI::Acquisition::Seismogram<ValueType> synthetic(receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P));
        
    waveconv_p_sum.assign(0);
    
   for (IndexType shotNumber = 0; shotNumber < sources.getNumShots(); shotNumber++) {
       
       /* --------------------------------------------------------------------------- */       
       /* Calculate wavefields and save them for every timestep in a #WavefieldRecord */
       /* --------------------------------------------------------------------------- */
       
       HOST_PRINT(comm, "\n================Start Forward====================\n");
       HOST_PRINT(comm, "Start time stepping for shot " << shotNumber + 1 << " of " << sources.getNumShots() << "\n"
                                                             << "Total Number of time steps: " << getNT << "\n");
       wavefields->reset();
       sources.init(config, ctx, dist, shotNumber);
       
       start_t = scai::common::Walltime::get();
       for (t = 0; t < getNT; t++) {

           solver.run(receivers, sources, model, *wavefields, derivatives, t, t + 1, config.get<ValueType>("DT"));

           // save wavefields in Dense Matrix
	   wavefieldrecord[t]->assign(*wavefields);

       }
       
       receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename") +".It" + std::to_string(iteration) +".shot" + std::to_string(shotNumber));
       
       end_t = scai::common::Walltime::get();
       HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");
       
       /* -------------------------------------------------------------------- */
       /* Calculate misfit and set seismograms of adjoint sources to residuals */
       /* -------------------------------------------------------------------- */
       
       std::string fieldSeisName("seismograms/rectangle.true");
       truedata.readFromFileRaw(fieldSeisName +".It0" + ".shot" + std::to_string(shotNumber) + ".p.mtx", adjoint.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P).getData().getRowDistributionPtr(), NULL);
       synthetic = receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P);
       
       synthetic -= truedata;
       misfitTemp.setValue(shotNumber, 0.5*synthetic.getData().l2Norm());  // misfit of one shot
       

       adjoint.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P) = synthetic;
       //  adjoint.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).writeToFileRaw("seismograms/adjoint.mtx");
       
       /* ------------------------------------------------------ */
       /*                Backward Modelling                      */
       /* ------------------------------------------------------ */

       HOST_PRINT(comm, "\n================Start Backward====================\n");
       HOST_PRINT(comm, "Start time stepping\n"
                               << "Total Number of time steps: " << getNT << "\n");
       
       wavefields->reset();
       waveconv_p = 0;
       waveconv_v = 0;
       
       start_t = scai::common::Walltime::get();

       for (t = getNT - 1; t >= 0; t--) {

           solver.run(receivers, adjoint, model, *wavefields, derivatives, t, t + 1, config.get<ValueType>("DT"));

           /* --------------------------------------- */
           /*             Convolution                 */
           /* --------------------------------------- */

	   tmp=wavefieldrecord[t]->getP();
           tmp *= wavefields->getP();
           waveconv_p += tmp;

//            wavefieldrecordvx.getColumn(tmp, t);
//            tmp *= wavefields->getVX();
//            waveconv_v += tmp;
//            wavefieldrecordvy.getColumn(tmp, t);
//            tmp *= wavefields->getVY();
//            waveconv_v += tmp;
        }
    
        waveconv_p_sum += waveconv_p;
//        waveconv_v_sum += waveconv_v;

        end_t = scai::common::Walltime::get();
        HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");

        //   receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).writeToFileRaw("seismograms/rec_adjoint.mtx");
        //
        //     waveconv_p.writeToFile(convname + ".shot_" + std::to_string(shotNumber)+".mtx");
        //=========================================================================================================

        } // end loop over shots
                
        dataMisfit.add(misfitTemp);
        /* Output jacobi */
        waveconv_p_sum.writeToFile(convname + "_p" +".It"+std::to_string(iteration)+ ".mtx");
//         waveconv_v_sum.writeToFile(convname + "_v" + ".mtx");

        /* ---------------------------------- */        
        /*       Calculate gradients          */
        /* ---------------------------------- */
        
        grad_bulk = model.getPWaveModulus();
        grad_bulk *= grad_bulk;
        grad_bulk.invert();
        grad_bulk *= waveconv_p_sum;
        grad_bulk *= -config.get<ValueType>("DT");

        grad_vp = 2 * grad_bulk;
        grad_vp *= model.getDensity();
        grad_vp *= model.getVelocityP();

//         grad_rho = model.getVelocityP();
//         grad_rho *= model.getVelocityP();
//         grad_rho *= grad_bulk;
//         grad_rho -= config.get<ValueType>("DT") * waveconv_v_sum;

        grad_vp.writeToFile(gradname + "_vp" + ".It"+std::to_string(iteration)+ ".mtx");
//         grad_rho.writeToFile(gradname + "_rho" + ".It"+std::to_string(iteration)+ ".mtx");
    

}

template class GradientCalculation<double>;
template class GradientCalculation<float>;




