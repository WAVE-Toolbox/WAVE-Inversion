#include <iostream>

#include "GradientCalculation.hpp"

template <typename ValueType>
void GradientCalculation<ValueType>::allocate(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr dist,  scai::hmemo::ContextPtr ctx)
{

    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    IndexType getNT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    wavefields = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefields->init(ctx, dist);

    ZeroLagXcorr = KITGPI::ZeroLagXcorr::Factory<ValueType>::Create(dimension, equationType);
    ZeroLagXcorr->init(ctx, dist);

    GradientPerShot = KITGPI::Gradient::Factory<ValueType>::Create(equationType);
    GradientPerShot->init(ctx, dist);
    /* ------------------------------------------------------- */
    /* Allocate wavefield record                               */
    /* ------------------------------------------------------- */

    for (IndexType i = 0; i < getNT; i++) {
        wavefieldPtr wavefieldsTemp(KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType));
        wavefieldsTemp->init(ctx, dist);

        wavefieldrecord.push_back(wavefieldsTemp);
    }
}

template <typename ValueType>
void GradientCalculation<ValueType>::calc(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Gradient::Gradient<ValueType> &gradient, KITGPI::Configuration::Configuration config, IndexType iteration, KITGPI::Misfit::Misfit<ValueType> &dataMisfit)
{

    double start_t, end_t; /* For timing */
    IndexType t = 0;
    IndexType getNT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */

    scai::dmemo::DistributionPtr dist = wavefields->getRefVX().getDistributionPtr();
    scai::dmemo::CommunicatorPtr comm = scai::dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                 // default context, set by environment variable SCAI_CONTEXT

    /* Should misfit be part of #GradientCalculation in a later stage? */
    scai::lama::DenseVector<ValueType> misfitTemp(sources.getNumShots(), 0, ctx);

    /* ------------------------------------------------------ */
    /* Allocate adjoint sources, true data and synthetic data */
    /* ------------------------------------------------------ */

    KITGPI::Acquisition::Receivers<ValueType> adjoint(config, ctx, dist);
    KITGPI::Acquisition::Seismogram<ValueType> truedata(adjoint.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P));
    KITGPI::Acquisition::Seismogram<ValueType> synthetic(receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P));

    //reset gradient because gradient is a sum of all gradientsPerShot gradients+=gradientPerShot
    gradient.reset();

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

        receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename") + ".It" + std::to_string(iteration) + ".shot" + std::to_string(shotNumber));

        end_t = scai::common::Walltime::get();
        HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");

        /* -------------------------------------------------------------------- */
        /* Calculate misfit and set seismograms of adjoint sources to residuals */
        /* -------------------------------------------------------------------- */

        std::string fieldSeisName(config.get<std::string>("FieldSeisName"));
        truedata.readFromFileRaw(fieldSeisName + ".It0" + ".shot" + std::to_string(shotNumber) + ".p.mtx", adjoint.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P).getData().getRowDistributionPtr(), NULL);
        synthetic = receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P);

        synthetic -= truedata;
        misfitTemp.setValue(shotNumber, 0.5 * synthetic.getData().l2Norm()); // misfit of one shot

        adjoint.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P) = synthetic;
        //  adjoint.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).writeToFileRaw("seismograms/adjoint.mtx");

        /* ------------------------------------------------------ */
        /*                Backward Modelling                      */
        /* ------------------------------------------------------ */

        HOST_PRINT(comm, "\n================Start Backward====================\n");
        HOST_PRINT(comm, "Start time stepping\n"
                             << "Total Number of time steps: " << getNT << "\n");

        wavefields->reset();
        ZeroLagXcorr->reset();

        start_t = scai::common::Walltime::get();

        for (t = getNT - 1; t >= 0; t--) {

            solver.run(receivers, adjoint, model, *wavefields, derivatives, t, t + 1, config.get<ValueType>("DT"));

            /* --------------------------------------- */
            /*             Convolution                 */
            /* --------------------------------------- */
            ZeroLagXcorr->update(*wavefieldrecord[t], *wavefields);
        }
        end_t = scai::common::Walltime::get();
        HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");

        /* ---------------------------------- */
        /*       Calculate gradients          */
        /* ---------------------------------- */
        GradientPerShot->estimateParameter(*ZeroLagXcorr, model, config.get<ValueType>("DT"));
        gradient += *GradientPerShot;

        //   receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).writeToFileRaw("seismograms/rec_adjoint.mtx");

    } // end loop over shots

    dataMisfit.add(misfitTemp);
}

template class GradientCalculation<double>;
template class GradientCalculation<float>;
