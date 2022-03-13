#include <iostream>
#include "GradientCalculation.hpp"

/*! \brief Allocation of wavefields, zero-lag crosscorrelation and gradients
 *
 \param config Configuration
 \param dist Distribution of the wave fields
 \param ctx Context
 \param workflow Workflow
 */
template <typename ValueType>
void KITGPI::GradientCalculation<ValueType>::allocate(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distInversion, scai::hmemo::ContextPtr ctx, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType numShotPerSuperShot)
{
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);   
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);

    scai::IndexType numRelaxationMechanisms = config.get<IndexType>("numRelaxationMechanisms");
    wavefields = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefields->init(ctx, dist, numRelaxationMechanisms);
    wavefieldsReflect = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefieldsReflect->init(ctx, dist, numRelaxationMechanisms);
    wavefieldsTemp = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefieldsTemp->init(ctx, distInversion, numRelaxationMechanisms);
    wavefieldsAdjointTemp = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefieldsAdjointTemp->init(ctx, distInversion, numRelaxationMechanisms);

    ZeroLagXcorr = KITGPI::ZeroLagXcorr::Factory<ValueType>::Create(dimension, equationType);
    ZeroLagXcorr->init(ctx, distInversion, workflow, config, numShotPerSuperShot);
    IndexType gradientKernel = config.getAndCatch("gradientKernel", 0); 
    IndexType decomposition = config.getAndCatch("decomposition", 0); 
    if ((gradientKernel == 2 || gradientKernel == 3) && decomposition == 0) {
        ZeroLagXcorrReflect = KITGPI::ZeroLagXcorr::Factory<ValueType>::Create(dimension, equationType);
        ZeroLagXcorrReflect->init(ctx, distInversion, workflow, config, numShotPerSuperShot);
    }
}

/*! \brief gather wavefields
 *
 \param wavefieldrecord Record of the wave fields
 \param wavefieldsInversion wavefields with a reduced size
 \param sourceFC frequency
 \param distInversion distribution of the wave fields
 \param t time in second
 */
template <typename ValueType>
void KITGPI::GradientCalculation<ValueType>::gatherWavefields(KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsInversion, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint, bool isReflect)
{
    if (!isReflect)
        ZeroLagXcorr->gatherWavefields(wavefieldsInversion, sourceFC, workflow, tStep, DT, isAdjoint);
    else
        ZeroLagXcorrReflect->gatherWavefields(wavefieldsInversion, sourceFC, workflow, tStep, DT, isAdjoint);
}

/*! \brief Initialization of the boundary conditions
 *
 *
 \param solver Forward solver
 \param derivatives Derivatives matrices
 \param receivers Receivers
 \param sources Sources 
 \param model Model for the finite-difference simulation
 \param gradientPerShot Gradient for simulations
 \param config Configuration
 \param dataMisfit Misfit
 */
template <typename ValueType>
void KITGPI::GradientCalculation<ValueType>::run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> sources, KITGPI::Acquisition::Receivers<ValueType> const &adjointSources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, std::vector<typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr> &wavefieldrecord, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, int shotNumber, int shotIndTrue, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Taper::Taper2D<ValueType> wavefieldTaper2D, std::vector<typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr> &wavefieldrecordReflect, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Preconditioning::EnergyPreconditioning<ValueType> energyPrecond, KITGPI::Preconditioning::EnergyPreconditioning<ValueType> energyPrecondReflect)
{
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
    ValueType DTinv = 1.0 / config.get<ValueType>("DT");
    double start_t_shot, end_t_shot; /* For timing */
    start_t_shot = common::Walltime::get();

    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower); 
    bool isSeismic = Common::checkEquationType<ValueType>(equationType);   
    scai::dmemo::DistributionPtr dist;
    scai::dmemo::CommunicatorPtr commShot;
    if (isSeismic) {
        if(equationType.compare("sh") == 0 || equationType.compare("viscosh") == 0){
            dist = wavefields->getRefVZ().getDistributionPtr();
        } else {
            dist = wavefields->getRefVX().getDistributionPtr();        
        }
        commShot = model.getDensity().getDistributionPtr()->getCommunicatorPtr(); // get communicator for shot domain
    } else {
        if(equationType.compare("tmem") == 0 || equationType.compare("viscotmem") == 0){
            dist = wavefields->getRefEZ().getDistributionPtr();
        } else {
            dist = wavefields->getRefEX().getDistributionPtr();        
        }
        commShot = model.getDielectricPermittivity().getDistributionPtr()->getCommunicatorPtr(); // get communicator for shot domain
    }
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                 // default context, set by environment variable SCAI_CONTEXT
    scai::dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    /* ------------------------------------------------------ */
    /*                Backward Modelling                      */
    /* ------------------------------------------------------ */
    IndexType gradientKernel = config.getAndCatch("gradientKernel", 0); 
    IndexType gradientDomain = config.getAndCatch("gradientDomain", 0);
    IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    IndexType decomposition = config.getAndCatch("decomposition", 0); 
    IndexType snapType = config.getAndCatch("snapType", 0);
    if (gradientKernel == 3) {
        IndexType numSwitch = gradientKernel - 2; 
        if ((workflow.iteration / numSwitch) % 2 == 0) {
            gradientKernel = 1;
        } else {
            gradientKernel = 2;
        }                
    } 
    if (decomposition != 0) {
        snapType = decomposition + 3;
    }  
    
    /* --------------------------------------- */
    /* Adjoint Wavefield record                */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> adjointSourcesReflect;   
    scai::dmemo::DistributionPtr distInversion = nullptr;   
    if (isSeismic) {
        if(equationType.compare("sh") == 0 || equationType.compare("viscosh") == 0){
            distInversion = wavefieldrecord[0]->getRefVZ().getDistributionPtr();
        } else {
            distInversion = wavefieldrecord[0]->getRefVX().getDistributionPtr();        
        }  
    } else {
        if(equationType.compare("tmem") == 0 || equationType.compare("viscotmem") == 0){
            distInversion = wavefieldrecord[0]->getRefEZ().getDistributionPtr();
        } else {
            distInversion = wavefieldrecord[0]->getRefEX().getDistributionPtr();        
        }  
    }          
    if (gradientKernel == 2 && decomposition == 0) { 
        adjointSourcesReflect.initWholeSpace(config, modelCoordinates, ctx, dist, receivers.getSeismogramTypes());
        ZeroLagXcorrReflect->prepareForInversion(gradientKernel, config);
    }     
    typename ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::SourceReceiverImplPtr SourceReceiverReflect(ForwardSolver::SourceReceiverImpl::Factory<ValueType>::Create(dimension, equationType, sources, adjointSourcesReflect, *wavefieldsReflect));
                
    wavefields->resetWavefields();
    ZeroLagXcorr->prepareForInversion(gradientKernel, config);
    bool isReflect = true;
    bool isAdjoint = true;
    
    lama::DenseVector<ValueType> compensation;
    if (config.getAndCatch("compensation", 0))
        compensation = model.getCompensation(config.get<ValueType>("DT"), 1);
    
    for (IndexType tStep = tStepEnd - 1; tStep > 0; tStep--) {
        *wavefieldsReflect = *wavefields;

        solver.run(receivers, adjointSources, model, *wavefields, derivatives, tStep);

        if (config.getAndCatch("compensation", 0))
            *wavefields *= compensation;
                
        if ((gradientKernel == 2 && decomposition == 0) || decomposition != 0) { 
            //calculate temporal derivative of wavefield
            *wavefieldsReflect -= *wavefields;
            *wavefieldsReflect *= -DTinv; // wavefieldsReflect will be gathered by adjointSourcesReflect
            if (gradientKernel == 2 && decomposition == 0) 
                SourceReceiverReflect->gatherSeismogram(tStep);
            if (decomposition != 0) 
                wavefields->decompose(decomposition, *wavefieldsReflect, derivatives);
        }
            
        if (((gradientKernel != 2 && decomposition == 0) || decomposition != 0) && tStep % workflow.skipDT == 0 && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep < tStepEnd / 2))) {
            wavefieldsAdjointTemp->applyTransform(wavefieldTaper2D.getAverageMatrix(), *wavefields);
            energyPrecond.intSquaredWavefields(*wavefieldsAdjointTemp, isAdjoint, config.get<ValueType>("DT"));
            if (gradientDomain == 0) { 
                /*  Cross correlation in the time domain   */
                //calculate temporal derivative of wavefield
                *wavefieldsTemp = *wavefieldrecord[floor(tStep / workflow.skipDT + 0.5)];
                *wavefieldsTemp -= *wavefieldrecord[floor(tStep / workflow.skipDT - 0.5)];
                *wavefieldsTemp *= DTinv;      
                *wavefieldsTemp *= workflow.skipDT; 
            
                /* please note that we exchange the position of the derivative and the forwardwavefield itself, which is different with the defination in ZeroLagXcorr function */
                ZeroLagXcorr->update(*wavefieldsTemp, *wavefieldrecord[floor(tStep / workflow.skipDT + 0.5)], *wavefieldsAdjointTemp, workflow);
                
                if (workflow.workflowStage == 0 && config.getAndCatch("snapType", 0) > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT")) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), config.get<ValueType>("DT")) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT"))) % Common::time2index(config.get<ValueType>("tincSnapshot"), config.get<ValueType>("DT")) == 0) {
                    if (gradientDomain == 0) {
                        ZeroLagXcorr->write(config.getAndCatch<std::string>("WavefieldFileName", "") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber), tStep, workflow);
                    }
                    wavefields->write(snapType, config.getAndCatch<std::string>("WavefieldFileName", "") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".receiver", tStep, derivatives, model, config.get<IndexType>("FileFormat"));
                }
            } else if (gradientDomain == 1 || gradientDomain == 2) {
                /* Cross correlation in the frequency domain */
                ZeroLagXcorr->gatherWavefields(*wavefieldsAdjointTemp, sources.getSourceFC(shotIndTrue), workflow, tStep, config.get<ValueType>("DT"), isAdjoint);
            } else if (gradientDomain == 3 && tStep < tStepEnd / 2) {
                this->gatherWavefields(*wavefieldsAdjointTemp, sources.getSourceFC(shotIndTrue), workflow, tStep, config.get<ValueType>("DT"), isAdjoint);
            }
        } else if (gradientKernel == 2 && decomposition == 0 && tStep % workflow.skipDT == 0 && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep < tStepEnd / 2))) {
            wavefieldsAdjointTemp->applyTransform(wavefieldTaper2D.getAverageMatrix(), *wavefields);
            energyPrecond.intSquaredWavefields(*wavefieldsAdjointTemp, isAdjoint, config.get<ValueType>("DT"));
            if (gradientDomain == 0) { 
                /*  Cross correlation in the time domain   */
                //calculate temporal derivative of wavefield
                *wavefieldsTemp = *wavefieldrecordReflect[floor(tStep / workflow.skipDT + 0.5)];
                *wavefieldsTemp -= *wavefieldrecordReflect[floor(tStep / workflow.skipDT - 0.5)];
                *wavefieldsTemp *= DTinv;      
                *wavefieldsTemp *= workflow.skipDT; 
            
                /* please note that we exchange the position of the derivative and the forwardwavefield itself, which is different with the defination in ZeroLagXcorr function */
                ZeroLagXcorrReflect->update(*wavefieldsTemp, *wavefieldrecordReflect[floor(tStep / workflow.skipDT + 0.5)], *wavefieldsAdjointTemp, workflow);
                
                if (workflow.workflowStage == 0 && config.getAndCatch("snapType", 0) > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT")) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), config.get<ValueType>("DT")) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT"))) % Common::time2index(config.get<ValueType>("tincSnapshot"), config.get<ValueType>("DT")) == 0) {
                    if (gradientDomain == 0) {
                        ZeroLagXcorrReflect->write(config.getAndCatch<std::string>("WavefieldFileName", "") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".sourceReflect", tStep, workflow);
                    }
                    wavefields->write(snapType, config.getAndCatch<std::string>("WavefieldFileName", "") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".receiver", tStep, derivatives, model, config.get<IndexType>("FileFormat"));
                }
            } else if (gradientDomain == 1 || gradientDomain == 2) {
                /* Cross correlation in the frequency domain */
                ZeroLagXcorrReflect->gatherWavefields(*wavefieldsAdjointTemp, sources.getSourceFC(shotIndTrue), workflow, tStep, config.get<ValueType>("DT"), isAdjoint);
            } else if (gradientDomain == 3 && tStep < tStepEnd / 2) {
                this->gatherWavefields(*wavefieldsAdjointTemp, sources.getSourceFC(shotIndTrue), workflow, tStep, config.get<ValueType>("DT"), isAdjoint, isReflect);
            }
        }
    }
    solver.resetCPML();

    // check wavefield for NaNs or infinite values
    if (commShot->any(!wavefields->isFinite(dist)) && commInterShot->getRank()==0){ // if any processor returns isfinite=false, write model and break
        model.write("model_crash", config.get<IndexType>("FileFormat"));
        COMMON_THROWEXCEPTION("Infinite or NaN value in adjoint wavefield, output model as model_crash.FILE_EXTENSION!");
    }

    /* ---------------------------------- */
    /*       Calculate gradients          */
    /* ---------------------------------- */
    if (gradientDomain != 0) {
        /* Cross correlation in the frequency domain */
        std::string filename = config.getAndCatch<std::string>("WavefieldFileName", "") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber);
        if (gradientKernel == 2 && decomposition == 0) {
            filename += ".sourceReflect";
        }
        ZeroLagXcorr->sumWavefields(commShot, filename, config.getAndCatch("snapType", 0), workflow, sources.getSourceFC(shotIndTrue), config.get<ValueType>("DT"), shotNumber);
    }
    ZeroLagXcorr->applyTransform(wavefieldTaper2D.getRecoverMatrix(), workflow);
    gradientPerShot.estimateParameter(*ZeroLagXcorr, model, config.get<ValueType>("DT"), workflow);
    ZeroLagXcorr->resetXcorr(workflow);
    
    /* Apply receiver Taper (if sourceTaperRadius=0 gradient will be multiplied by 1) */
    SourceTaper.init(dist, ctx, sources, config, modelCoordinates, config.get<IndexType>("sourceTaperRadius"));
    SourceTaper.apply(gradientPerShot);    
    /* Apply receiver Taper (if receiverTaperRadius=0 gradient will be multiplied by 1) */
    ReceiverTaper.init(dist, ctx, receivers, config, modelCoordinates, config.get<IndexType>("receiverTaperRadius"));
    ReceiverTaper.apply(gradientPerShot);
    sourceReceiverTaper.init(dist, ctx, sources, receivers, config, modelCoordinates);
    sourceReceiverTaper.apply(gradientPerShot);

    /* Apply energy preconditioning per shot */
    energyPrecond.applyTransform(wavefieldTaper2D.getRecoverMatrix());
    energyPrecond.apply(gradientPerShot, shotNumber, config.get<IndexType>("FileFormat"));
    gradientPerShot.applyMedianFilter(commAll, config);   
    
    scai::lama::DenseVector<ValueType> mask; //mask to restore vacuum
    if (isSeismic) {
        if(equationType.compare("sh") == 0 || equationType.compare("viscosh") == 0){
            mask = model.getVelocityS();  
        } else {
            mask = model.getVelocityP();      
        }  
    } else {    
        mask = model.getDielectricPermittivity();
        mask /= model.getDielectricPermittivityVacuum();  // calculate the relative dielectricPermittivity    
        mask -= 1;
    }
    mask.unaryOp(mask, common::UnaryOp::SIGN);
    mask.unaryOp(mask, common::UnaryOp::ABS); 
    gradientPerShot *= mask;

    gradientPerShot.normalize();
    
    end_t_shot = common::Walltime::get();
    HOST_PRINT(commShot, "Shot number " << shotNumber << ": Finish gradient calculation in " << end_t_shot - start_t_shot << " sec.\n");
            
    if (gradientKernel == 2 && decomposition == 0) {
        typename KITGPI::Gradient::Gradient<ValueType>::GradientPtr testgradient(KITGPI::Gradient::Factory<ValueType>::Create(equationType));
        *testgradient = gradientPerShot;
        scai::lama::DenseVector<ValueType> reflectivity;
        reflectivity = model.getReflectivity();
        dataMisfit.calcReflectSources(adjointSourcesReflect, reflectivity);
        wavefields->resetWavefields();
    
        for (IndexType tStep = tStepEnd - 1; tStep > 0; tStep--) {
            
            solver.run(receivers, adjointSourcesReflect, model, *wavefields, derivatives, tStep);

            if (config.getAndCatch("compensation", 0))
                *wavefields *= compensation;
            
            /* --------------------------------------- */
            /*  Cross correlation in the time domain   */
            /* --------------------------------------- */
            if (tStep % workflow.skipDT == 0 && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep < tStepEnd / 2))) {
                wavefieldsAdjointTemp->applyTransform(wavefieldTaper2D.getAverageMatrix(), *wavefields);
                energyPrecondReflect.intSquaredWavefields(*wavefieldsAdjointTemp, isAdjoint, config.get<ValueType>("DT"));
                if (gradientDomain == 0) { 
                    /*  Cross correlation in the time domain   */
                    //calculate temporal derivative of wavefield
                    *wavefieldsTemp = *wavefieldrecord[floor(tStep / workflow.skipDT + 0.5)];
                    *wavefieldsTemp -= *wavefieldrecord[floor(tStep / workflow.skipDT - 0.5)];
                    *wavefieldsTemp *= DTinv;      
                    *wavefieldsTemp *= workflow.skipDT; 
                
                    /* please note that we exchange the position of the derivative and the forwardwavefield itself, which is different with the defination in ZeroLagXcorr function */
                    ZeroLagXcorr->update(*wavefieldsTemp, *wavefieldrecord[floor(tStep / workflow.skipDT + 0.5)], *wavefieldsAdjointTemp, workflow);
                    
                    if (workflow.workflowStage == 0 && config.getAndCatch("snapType", 0) > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT")) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), config.get<ValueType>("DT")) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT"))) % Common::time2index(config.get<ValueType>("tincSnapshot"), config.get<ValueType>("DT")) == 0) {
                        if (gradientDomain == 0) {
                            ZeroLagXcorr->write(config.getAndCatch<std::string>("WavefieldFileName", "") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".receiverReflect", tStep, workflow);
                        }
                        wavefields->write(config.getAndCatch("snapType", 0), config.getAndCatch<std::string>("WavefieldFileName", "") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".receiverReflect", tStep, derivatives, model, config.get<IndexType>("FileFormat"));
                    }
                } else if (gradientDomain == 1 || gradientDomain == 2) {
                    /* Cross correlation in the frequency domain */
                    ZeroLagXcorr->gatherWavefields(*wavefieldsAdjointTemp, sources.getSourceFC(shotIndTrue), workflow, tStep, config.get<ValueType>("DT"), isAdjoint);
                } else if (gradientDomain == 3 && tStep < tStepEnd / 2) {
                    this->gatherWavefields(*wavefieldsAdjointTemp, sources.getSourceFC(shotIndTrue), workflow, tStep, config.get<ValueType>("DT"), isAdjoint);
                }
            }
        }       
        solver.resetCPML();
        
        /* ---------------------------------- */
        /*       Calculate gradients          */
        /* ---------------------------------- */ 
        if (gradientDomain != 0) {
            /* Cross correlation in the frequency domain */
            ZeroLagXcorrReflect->sumWavefields(commShot, config.getAndCatch<std::string>("WavefieldFileName", "") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".receiverReflect", config.getAndCatch("snapType", 0), workflow, sources.getSourceFC(shotIndTrue), config.get<ValueType>("DT"), shotNumber);
        }  
        ZeroLagXcorrReflect->applyTransform(wavefieldTaper2D.getRecoverMatrix(), workflow); 
        gradientPerShot.estimateParameter(*ZeroLagXcorrReflect, model, config.get<ValueType>("DT"), workflow);
        ZeroLagXcorrReflect->resetXcorr(workflow);
    
        SourceTaper.apply(gradientPerShot);
        ReceiverTaper.apply(gradientPerShot);
        sourceReceiverTaper.apply(gradientPerShot);

        /* Apply energy preconditioning per shot */
        energyPrecondReflect.applyTransform(wavefieldTaper2D.getRecoverMatrix());
        energyPrecondReflect.apply(gradientPerShot, shotNumber, config.get<IndexType>("FileFormat"));
        gradientPerShot.applyMedianFilter(commAll, config); 
        gradientPerShot *= mask;
        gradientPerShot.normalize(); 
        gradientPerShot += *testgradient;
        
        end_t_shot = common::Walltime::get();
        HOST_PRINT(commShot, "Shot number " << shotNumber << ": Finish tomographic gradient calculation in " << end_t_shot - start_t_shot << " sec.\n");
    }
    
    if (config.get<IndexType>("writeGradientPerShot"))
        gradientPerShot.write(config.get<std::string>("GradientFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber), config.get<IndexType>("FileFormat"), workflow);
}

template class KITGPI::GradientCalculation<double>;
template class KITGPI::GradientCalculation<float>;
