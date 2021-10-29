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
void KITGPI::GradientCalculation<ValueType>::run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> const &adjointSources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Gradient::Gradient<ValueType> &gradient, std::vector<typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr> &wavefieldrecord, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, int shotNumber, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Taper::Taper2D<ValueType> wavefieldTaper2D, std::vector<typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr> &wavefieldrecordReflect, KITGPI::Misfit::Misfit<ValueType> &dataMisfit)
{
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
    IndexType dtinversion = config.get<IndexType>("DTInversion");    
    ValueType DTinv = 1 / config.get<ValueType>("DT");

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
    IndexType gradientType = config.getAndCatch("gradientType", 0); 
    IndexType decomposeType = config.getAndCatch("decomposeType", 0); 
    IndexType snapType = config.getAndCatch("snapType", 0);
    if (gradientType > 2) {
        IndexType numSwitch = gradientType - 2; 
        if ((workflow.iteration / numSwitch) % 2 == 0) {
            gradientType = 1;
        } else {
            gradientType = 2;
        }                
    } 
    if (decomposeType != 0) {
        snapType = decomposeType + 3;
    }  
    
    energyPrecond.init(dist, config);
    wavefields->resetWavefields();
    ZeroLagXcorr->setDecomposeType(decomposeType);
    ZeroLagXcorr->setGradientType(gradientType);
    ZeroLagXcorr->init(ctx, dist, workflow);
    
    /* --------------------------------------- */
    /* Adjoint Wavefield record                */
    /* --------------------------------------- */
    std::vector<wavefieldPtr> wavefieldrecordAdjointReflect;  
    Acquisition::Receivers<ValueType> adjointSourcesReflect;                
    if (gradientType == 2 && decomposeType == 0) { 
        scai::dmemo::DistributionPtr distInversion;
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
        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
            if (tStep % dtinversion == 0) {
                wavefieldPtr wavefieldsInversion = Wavefields::Factory<ValueType>::Create(dimension, equationType);
                wavefieldsInversion->init(ctx, distInversion);
                wavefieldrecordAdjointReflect.push_back(wavefieldsInversion);
            }
        }
        adjointSourcesReflect.initWholeSpace(config, modelCoordinates, ctx, dist, receivers.getSeismogramTypes());
    }     
    typename ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::SourceReceiverImplPtr SourceReceiverReflect(ForwardSolver::SourceReceiverImpl::Factory<ValueType>::Create(dimension, equationType, sources, adjointSourcesReflect, *wavefieldsTemp));
    
    for (IndexType tStep = tStepEnd - 1; tStep > 0; tStep--) {
        *wavefieldsTemp = *wavefields;

        solver.run(receivers, adjointSources, model, *wavefields, derivatives, tStep);

        if ((gradientType == 2 && decomposeType == 0) || decomposeType != 0) { 
            //calculate temporal derivative of wavefield
            *wavefieldsTemp -= *wavefields;
            *wavefieldsTemp *= -DTinv; // wavefieldsTemp will be gathered by adjointSourcesReflect
            if (gradientType == 2 && decomposeType == 0) 
                SourceReceiverReflect->gatherSeismogram(tStep);
            if (decomposeType != 0) 
                wavefields->decompose(decomposeType, *wavefieldsTemp, derivatives);
            /* --------------------------------------- */
            /*             Convolution                 */
            /* --------------------------------------- */
            if (gradientType == 2 && decomposeType == 0 && tStep % dtinversion == 0) {
                // save wavefields in std::vector
                *wavefieldrecordAdjointReflect[floor(tStep / dtinversion + 0.5)] = wavefieldTaper2D.applyWavefieldAverage(wavefields);
                *wavefieldsTemp = wavefieldTaper2D.applyWavefieldRecover(wavefieldrecordReflect[floor(tStep / dtinversion + 0.5)]);
                energyPrecond.intSquaredWavefields(*wavefieldsTemp, *wavefields, config.get<ValueType>("DT"));
                //calculate temporal derivative of wavefield
                *wavefieldsTemp -= wavefieldTaper2D.applyWavefieldRecover(wavefieldrecordReflect[floor(tStep / dtinversion - 0.5)]);
                *wavefieldsTemp *= DTinv;      
                *wavefieldsTemp *= dtinversion; 
            
                /* please note that we exchange the position of the derivative and the forwardwavefield itself, which is different with the defination in ZeroLagXcorr function */
                ZeroLagXcorr->update(*wavefieldsTemp, wavefieldTaper2D.applyWavefieldRecover(wavefieldrecordReflect[floor(tStep / dtinversion + 0.5)]), *wavefields, workflow);
                
                if (workflow.workflowStage == 0 && config.getAndCatch("snapType", 0) > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT")) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), config.get<ValueType>("DT")) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT"))) % Common::time2index(config.get<ValueType>("tincSnapshot"), config.get<ValueType>("DT")) == 0) {
                    ZeroLagXcorr->write(config.get<std::string>("WavefieldFileName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".sourceReflect", tStep, workflow);
                    wavefields->write(snapType, config.get<std::string>("WavefieldFileName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".receiver", tStep, derivatives, model, config.get<IndexType>("FileFormat"));
                }
            }
        }        
        if (((gradientType != 2 && decomposeType == 0) || decomposeType != 0) && tStep % dtinversion == 0) {
            /* --------------------------------------- */
            /*             Convolution                 */
            /* --------------------------------------- */
            *wavefieldsTemp = wavefieldTaper2D.applyWavefieldRecover(wavefieldrecord[floor(tStep / dtinversion + 0.5)]);
            energyPrecond.intSquaredWavefields(*wavefieldsTemp, *wavefields, config.get<ValueType>("DT"));
            //calculate temporal derivative of wavefield
            *wavefieldsTemp -= wavefieldTaper2D.applyWavefieldRecover(wavefieldrecord[floor(tStep / dtinversion - 0.5)]);
            *wavefieldsTemp *= DTinv;      
            *wavefieldsTemp *= dtinversion; 
           
            /* please note that we exchange the position of the derivative and the forwardwavefield itself, which is different with the defination in ZeroLagXcorr function */
            ZeroLagXcorr->update(*wavefieldsTemp, wavefieldTaper2D.applyWavefieldRecover(wavefieldrecord[floor(tStep / dtinversion + 0.5)]), *wavefields, workflow);
            
            if (workflow.workflowStage == 0 && config.getAndCatch("snapType", 0) > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT")) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), config.get<ValueType>("DT")) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT"))) % Common::time2index(config.get<ValueType>("tincSnapshot"), config.get<ValueType>("DT")) == 0) {
                ZeroLagXcorr->write(config.get<std::string>("WavefieldFileName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + "", tStep, workflow);
                wavefields->write(snapType, config.get<std::string>("WavefieldFileName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".receiver", tStep, derivatives, model, config.get<IndexType>("FileFormat"));
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
    gradient.estimateParameter(*ZeroLagXcorr, model, config.get<ValueType>("DT"), workflow);
    SourceTaper.init(dist, ctx, sources, config, modelCoordinates, config.get<IndexType>("sourceTaperRadius"));
    SourceTaper.apply(gradient);
    /* Apply receiver Taper (if ReceiverTaperRadius=0 gradient will be multplied by 1) */
    ReceiverTaper.init(dist, ctx, receivers, config, modelCoordinates, config.get<IndexType>("receiverTaperRadius"));
    ReceiverTaper.apply(gradient);
    sourceReceiverTaper.init(dist, ctx, sources, receivers, config, modelCoordinates);
    sourceReceiverTaper.apply(gradient);

    /* Apply energy preconditioning per shot */
    energyPrecond.apply(gradient, shotNumber, config.get<IndexType>("FileFormat"));
    gradient.applyMedianFilter(config); 
    
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
    gradient *= mask;
//     IO::writeVector(mask, "gradients/mask", 1);

    gradient.normalize();
            
    if (gradientType == 2 && decomposeType == 0) {
        typename KITGPI::Gradient::Gradient<ValueType>::GradientPtr testgradient(KITGPI::Gradient::Factory<ValueType>::Create(equationType));
        *testgradient = gradient;
        scai::lama::DenseVector<ValueType> reflectivity;
        reflectivity = model.getReflectivity();
        dataMisfit.calcReflectSources(adjointSourcesReflect, reflectivity);
        energyPrecond.init(dist, config);
        wavefields->resetWavefields();
        ZeroLagXcorr->resetXcorr(workflow);
        for (IndexType tStep = tStepEnd - 1; tStep > 0; tStep--) {
            
            solver.run(receivers, adjointSourcesReflect, model, *wavefields, derivatives, tStep);

            /* --------------------------------------- */
            /*             Convolution                 */
            /* --------------------------------------- */
            if (tStep % dtinversion == 0) {
                *wavefieldsTemp = wavefieldTaper2D.applyWavefieldRecover(wavefieldrecord[floor(tStep / dtinversion + 0.5)]);
                energyPrecond.intSquaredWavefields(*wavefieldsTemp, *wavefields, config.get<ValueType>("DT"));
                //calculate temporal derivative of wavefield
                *wavefieldsTemp -= wavefieldTaper2D.applyWavefieldRecover(wavefieldrecord[floor(tStep / dtinversion - 0.5)]);
                *wavefieldsTemp *= DTinv;      
                *wavefieldsTemp *= dtinversion; 
            
                /* please note that we exchange the position of the derivative and the forwardwavefield itself, which is different with the defination in ZeroLagXcorr function */
                ZeroLagXcorr->update(*wavefieldsTemp, wavefieldTaper2D.applyWavefieldRecover(wavefieldrecord[floor(tStep / dtinversion + 0.5)]), *wavefields, workflow);
                
                if (workflow.workflowStage == 0 && config.getAndCatch("snapType", 0) > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT")) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), config.get<ValueType>("DT")) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT"))) % Common::time2index(config.get<ValueType>("tincSnapshot"), config.get<ValueType>("DT")) == 0) {
                    ZeroLagXcorr->write(config.get<std::string>("WavefieldFileName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".receiverReflect", tStep, workflow);
                    wavefields->write(config.getAndCatch("snapType", 0), config.get<std::string>("WavefieldFileName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".receiverReflect", tStep, derivatives, model, config.get<IndexType>("FileFormat"));
                }
            }
        }       
        solver.resetCPML();
        
        /* ---------------------------------- */
        /*       Calculate gradients          */
        /* ---------------------------------- */    
        gradient.estimateParameter(*ZeroLagXcorr, model, config.get<ValueType>("DT"), workflow);
        SourceTaper.apply(gradient);
        ReceiverTaper.apply(gradient);
        sourceReceiverTaper.apply(gradient);

        /* Apply energy preconditioning per shot */
        energyPrecond.apply(gradient, shotNumber, config.get<IndexType>("FileFormat"));
        gradient.applyMedianFilter(config); 
        gradient *= mask;
        gradient.normalize(); 
        gradient += *testgradient;
    }
    
    if (config.get<IndexType>("writeGradientPerShot"))
        gradient.write(config.get<std::string>("GradientFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber), config.get<IndexType>("FileFormat"), workflow);
}

template class KITGPI::GradientCalculation<double>;
template class KITGPI::GradientCalculation<float>;
