#include <iostream>

#include "GradientCalculation.hpp"

/*! \brief Allocation of wavefields, zero-lag crosscorrelation and gradients
 *
 *
 \param config Configuration
 \param dist Distribution of the wave fields
 \param ctx Context
 */
template <typename ValueType>
void GradientCalculation<ValueType>::allocate(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr dist,  scai::hmemo::ContextPtr ctx)
{

    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");

    wavefields = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefields->init(ctx, dist);

    ZeroLagXcorr = KITGPI::ZeroLagXcorr::Factory<ValueType>::Create(dimension, equationType);
    ZeroLagXcorr->init(config,ctx, dist);

}

/*! \brief Initialitation of the boundary conditions
 *
 *
 \param solver Forward solver
 \param derivatives Derivatives matrices
 \param receivers Receivers
 \param sources Sources 
 \param model Model for the finite-difference simulation
 \param gradient Gradient for simulations
 \param config Configuration
 \param iteration Current iteration count
 \param dataMisfit Misfit
 */
template <typename ValueType>
void GradientCalculation<ValueType>::run(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> const &adjointSources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Gradient::Gradient<ValueType> &gradient, std::vector<typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr> &wavefieldrecord, KITGPI::Configuration::Configuration config, IndexType iteration, int shotNumber)
{

    IndexType t = 0;
    IndexType tEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */

    scai::dmemo::DistributionPtr dist = wavefields->getRefVX().getDistributionPtr();
    scai::dmemo::CommunicatorPtr comm = scai::dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                 // default context, set by environment variable SCAI_CONTEXT

    /* ------------------------------------------------------ */
    /*                Backward Modelling                      */
    /* ------------------------------------------------------ */

    HOST_PRINT(comm, "\n-------------- Start Backward -------------------\n");

    wavefields->resetWavefields();
    ZeroLagXcorr->resetXcorr();


    for (t = tEnd - 1; t >= 0; t--) {

        solver.run(receivers, adjointSources, model, *wavefields, derivatives, t, t + 1, config.get<ValueType>("DT"));

        /* --------------------------------------- */
        /*             Convolution                 */
        /* --------------------------------------- */
        ZeroLagXcorr->update(*wavefieldrecord[t], *wavefields);
    }


    /* ---------------------------------- */
    /*       Calculate gradients          */
    /* ---------------------------------- */
    HOST_PRINT(comm, "\nCalculate Gradient\n");
    gradient.estimateParameter(*ZeroLagXcorr, model, config.get<ValueType>("DT"));
    SourceTaper.init(dist,ctx,sources,config,config.get<IndexType>("SourceTaperRadius"));
    SourceTaper.apply(gradient);

    if(config.get<IndexType>("WriteGradientPerShot"))
       gradient.getVelocityP().writeToFile(config.get<std::string>("GradientFilename")  + ".It_" + std::to_string(iteration + 1) + ".Shot_" + std::to_string(shotNumber) + ".vp" + ".mtx");

//     receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).writeToFileRaw("seismograms/rec_adjoint.mtx");
    
}

template class GradientCalculation<double>;
template class GradientCalculation<float>;
