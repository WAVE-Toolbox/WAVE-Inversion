#include <iostream>

#include "GradientCalculation.hpp"

using scai::IndexType;

/*! \brief Allocation of wavefieldsEM, zero-lag crosscorrelation and gradients
 *
 *
 \param configEM Configuration
 \param distEM Distribution of the wave fields
 \param ctx Context
 \param workflowEM Workflow
 */
template <typename ValueType>
void KITGPI::GradientCalculationEM<ValueType>::allocate(KITGPI::Configuration::Configuration configEM, scai::dmemo::DistributionPtr distEM, scai::hmemo::ContextPtr ctx, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    std::string dimensionEM = configEM.get<std::string>("dimension");
    std::string equationTypeEM = configEM.get<std::string>("equationType");
    std::transform(dimensionEM.begin(), dimensionEM.end(), dimensionEM.begin(), ::tolower);   
    std::transform(equationTypeEM.begin(), equationTypeEM.end(), equationTypeEM.begin(), ::tolower);   

    wavefieldsEM = KITGPI::Wavefields::FactoryEM<ValueType>::Create(dimensionEM, equationTypeEM);
    wavefieldsEM->init(ctx, distEM);
    wavefieldsTempEM = KITGPI::Wavefields::FactoryEM<ValueType>::Create(dimensionEM, equationTypeEM);
    wavefieldsTempEM->init(ctx, distEM);

    ZeroLagXcorr = KITGPI::ZeroLagXcorr::FactoryEM<ValueType>::Create(dimensionEM, equationTypeEM);
    ZeroLagXcorr->init(ctx, distEM, workflowEM);
}

/*! \brief Initialitation of the boundary conditions
 *
 *
 \param solverEM Forward solverEM
 \param derivativesEM Derivatives matrices
 \param ReceiversEM Receivers
 \param sourcesEM Sources 
 \param modelEM Model for the finite-difference simulation
 \param gradientEM Gradient for simulations
 \param configEM Configuration
 \param dataMisfitEM Misfit
 */
template <typename ValueType>
void KITGPI::GradientCalculationEM<ValueType>::run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolverEM<ValueType> &solverEM, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivativesEM, KITGPI::Acquisition::ReceiversEM<ValueType> &ReceiversEM, KITGPI::Acquisition::SourcesEM<ValueType> &sourcesEM, KITGPI::Acquisition::ReceiversEM<ValueType> &adjointSourcesEM, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Gradient::GradientEM<ValueType> &gradientEM, std::vector<typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr> &wavefieldrecordEM, KITGPI::Configuration::Configuration configEM, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesEM, int shotNumber, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM, KITGPI::Taper::Taper2D<ValueType> &gradientTaper2DPerShotEM)
{
    IndexType t = 0;
    IndexType tEnd = static_cast<IndexType>((configEM.get<ValueType>("T") / configEM.get<ValueType>("DT")) + 0.5);

    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    std::string equationTypeEM = configEM.get<std::string>("equationType");
    std::transform(equationTypeEM.begin(), equationTypeEM.end(), equationTypeEM.begin(), ::tolower);  
    scai::dmemo::DistributionPtr distEM;
    if(equationTypeEM.compare("tmem") == 0 || equationTypeEM.compare("viscotmem") == 0){
        distEM = wavefieldsEM->getRefEZ().getDistributionPtr();
    } else {
        distEM = wavefieldsEM->getRefEX().getDistributionPtr();        
    }
    scai::dmemo::CommunicatorPtr commShot = modelEM.getDielectricPermittivityEM().getDistributionPtr()->getCommunicatorPtr(); // get communicator for shot domain
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                 // default context, set by environment variable SCAI_CONTEXT
    scai::dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    /* ------------------------------------------------------ */
    /*                Backward Modelling                      */
    /* ------------------------------------------------------ */

    wavefieldsEM->resetWavefields();
    ZeroLagXcorr->resetXcorr(workflowEM);

    IndexType dtinversionEM = configEM.get<IndexType>("DTInversion");
    
    ValueType DTinv = 1 / configEM.get<ValueType>("DT");
    
    /* --------------------------------------- */
    /* Adjoint Wavefield record                        */
    /* --------------------------------------- */
    
//     std::cout << "Shot " << shotInd + 1 << ": Start adjoint solverEM\n" << std::endl;
    for (t = tEnd - 1; t >= 0; t--) {
        
        *wavefieldsTempEM = *wavefieldsEM;
        
        solverEM.run(ReceiversEM, adjointSourcesEM, modelEM, *wavefieldsEM, derivativesEM, t);

        /* --------------------------------------- */
        /*             Convolution                 */
        /* --------------------------------------- */
        if (t % dtinversionEM == 0) {
            if (t != 0) {
                *wavefieldsTempEM = gradientTaper2DPerShotEM.applyWavefieldRecover(wavefieldrecordEM[floor(t / dtinversionEM + 0.5)]);
                //calculate temporal derivative of wavefield
                *wavefieldsTempEM -= gradientTaper2DPerShotEM.applyWavefieldRecover(wavefieldrecordEM[floor(t / dtinversionEM - 0.5)]);
                *wavefieldsTempEM *= DTinv;      
            } else {
                *wavefieldsTempEM = gradientTaper2DPerShotEM.applyWavefieldRecover(wavefieldrecordEM[floor(t / dtinversionEM + 1.5)]);
                //calculate temporal derivative of wavefield
                *wavefieldsTempEM -= gradientTaper2DPerShotEM.applyWavefieldRecover(wavefieldrecordEM[floor(t / dtinversionEM + 0.5)]);
                *wavefieldsTempEM *= DTinv;      
            }
            *wavefieldsTempEM *= dtinversionEM; 
           
            /* please note that we exchange the position of the derivative and the forwardwavefield itself, which is different with the defination in ZeroLagXcorr function */
            ZeroLagXcorr->update(*wavefieldsTempEM, gradientTaper2DPerShotEM.applyWavefieldRecover(wavefieldrecordEM[floor(t / dtinversionEM + 0.5)]), *wavefieldsEM, workflowEM);
        }
    }

    // check wavefield for NaNs or infinite values
    if (commShot->any(!wavefieldsEM->isFinite(distEM)) && commInterShot->getRank()==0){ // if any processor returns isfinite=false, write model and break
        modelEM.write("model_crash", configEM.get<IndexType>("FileFormat"));
        COMMON_THROWEXCEPTION("Infinite or NaN value in adjoint EM wavefield, output model as model_crash.FILE_EXTENSION!");
    }
    
    /* ---------------------------------- */
    /*       Calculate gradients          */
    /* ---------------------------------- */
    
    gradientEM.estimateParameter(*ZeroLagXcorr, modelEM, configEM.get<ValueType>("DT"), workflowEM);
    SourceTaper.init(distEM, ctx, sourcesEM, configEM, modelCoordinatesEM, configEM.get<IndexType>("sourceTaperRadius"));
    SourceTaper.apply(gradientEM);

    if (configEM.get<IndexType>("writeGradientPerShot"))
        gradientEM.write(configEM.get<std::string>("GradientFilename") + ".stage_" + std::to_string(workflowEM.workflowStage + 1) + ".It_" + std::to_string(workflowEM.iteration + 1) + ".shot_" + std::to_string(shotNumber), configEM.get<IndexType>("FileFormat"), workflowEM);
}

template class KITGPI::GradientCalculationEM<double>;
template class KITGPI::GradientCalculationEM<float>;
