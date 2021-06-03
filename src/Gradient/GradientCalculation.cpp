#include <iostream>

#include "GradientCalculation.hpp"

using scai::IndexType;

/*! \brief Allocation of wavefields, zero-lag crosscorrelation and gradients
 *
 *
 \param config Configuration
 \param dist Distribution of the wave fields
 \param ctx Context
 \param workflow Workflow
 */
template <typename ValueType>
void KITGPI::GradientCalculation<ValueType>::allocate(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);   
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);

    wavefields = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefields->init(ctx, dist);
    wavefieldsTemp = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefieldsTemp->init(ctx, dist);

    ZeroLagXcorr = KITGPI::ZeroLagXcorr::Factory<ValueType>::Create(dimension, equationType);
    ZeroLagXcorr->init(ctx, dist, workflow);
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
 \param dataMisfit Misfit
 */
template <typename ValueType>
void KITGPI::GradientCalculation<ValueType>::run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> const &adjointSources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Gradient::Gradient<ValueType> &gradient, std::vector<typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr> &wavefieldrecord, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, int shotNumber, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Taper::Taper2D<ValueType> wavefieldTaper2D)
{
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    std::string equationType = config.get<std::string>("equationType");
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);  
    scai::dmemo::DistributionPtr dist;
    if(equationType.compare("sh") == 0){
        dist = wavefields->getRefVZ().getDistributionPtr();
    } else {
        dist = wavefields->getRefVX().getDistributionPtr();        
    }
    scai::dmemo::CommunicatorPtr commShot = model.getDensity().getDistributionPtr()->getCommunicatorPtr(); // get communicator for shot domain
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                 // default context, set by environment variable SCAI_CONTEXT
    scai::dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

 /* ------------------------------------------------------ */
    /*                Backward Modelling                      */
    /* ------------------------------------------------------ */

    wavefields->resetWavefields();
    ZeroLagXcorr->resetXcorr(workflow);

    IndexType dtinversion = config.get<IndexType>("DTInversion");
    ValueType DTinv = 1 / config.get<ValueType>("DT");

//     if (shotNumber == 0) {
//         lama::DenseVector<ValueType> temp1;
//         IndexType tStep = 1000;
//         *wavefieldsTemp = *wavefieldrecord[floor(tStep / dtinversion + 0.5)];
//         temp1 = wavefieldsTemp->getRefSxx();
//         std::cout << "\n GradientCalculation run wavefieldrecord[" << floor(tStep / dtinversion + 0.5) << "][6150]: " << temp1[6150] << "\n"<< std::endl;
//     }
    for (IndexType tStep = tStepEnd - 1; tStep > 0; tStep--) {

        solver.run(receivers, adjointSources, model, *wavefields, derivatives, tStep);

        /* --------------------------------------- */
        /*             Convolution                 */
        /* --------------------------------------- */
        if (tStep % dtinversion == 0) {
            *wavefieldsTemp = wavefieldTaper2D.applyWavefieldRecover(wavefieldrecord[floor(tStep / dtinversion + 0.5)]);
            //calculate temporal derivative of wavefield
            *wavefieldsTemp -= wavefieldTaper2D.applyWavefieldRecover(wavefieldrecord[floor(tStep / dtinversion - 0.5)]);
            *wavefieldsTemp *= DTinv;
            *wavefieldsTemp *= dtinversion; 
            
            ZeroLagXcorr->update(*wavefieldsTemp, wavefieldTaper2D.applyWavefieldRecover(wavefieldrecord[floor(tStep / dtinversion + 0.5)]), *wavefields, workflow);
        }
    }

    // check wavefield for NaNs or infinite values
    if (commShot->any(!wavefields->isFinite(dist)) && commInterShot->getRank()==0){ // if any processor returns isfinite=false, write model and break
        model.write("model_crash", config.get<IndexType>("FileFormat"));
        COMMON_THROWEXCEPTION("Infinite or NaN value in adjoint velocity wavefield, output model as model_crash.FILE_EXTENSION!");
    }

    /* ---------------------------------- */
    /*       Calculate gradients          */
    /* ---------------------------------- */
    gradient.estimateParameter(*ZeroLagXcorr, model, config.get<ValueType>("DT"), workflow);
    SourceTaper.init(dist, ctx, sources, config, modelCoordinates, config.get<IndexType>("sourceTaperRadius"));
    SourceTaper.apply(gradient);

    if (config.get<IndexType>("writeGradientPerShot"))
        gradient.write(config.get<std::string>("GradientFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber), config.get<IndexType>("FileFormat"), workflow);
}

template class KITGPI::GradientCalculation<double>;
template class KITGPI::GradientCalculation<float>;
