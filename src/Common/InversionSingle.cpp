#include "InversionSingle.hpp"

/*! \brief Constructor of the single inversion function
 \param config Configuration
 */
template <typename ValueType>
KITGPI::InversionSingle<ValueType>::InversionSingle(KITGPI::Configuration::Configuration config, IndexType inversionType)
{
    if (inversionType != 0) {
        dimension = config.get<std::string>("dimension");
        equationType = config.get<std::string>("equationType");
        misfitType = config.get<std::string>("misfitType");
        std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);   
        std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);
        std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::tolower);
        
        isSeismic = Common::checkEquationType<ValueType>(equationType);         
        useStreamConfig = config.getAndCatch("useStreamConfig", false);        
        if (useStreamConfig) {
            configBig.readFromFile(config.get<std::string>("streamConfigFilename"), true);
        }
        
        tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
        multiMisfitType = config.getAndCatch("multiMisfitType", misfitType);
        gradname = config.get<std::string>("gradientFilename");
        logFilename = config.get<std::string>("logFilename");
        steplengthInit = config.get<ValueType>("steplengthInit");
        optimizationType = config.get<std::string>("optimizationType");
        numRelaxationMechanisms = config.get<IndexType>("numRelaxationMechanisms");
        useSourceEncode = config.getAndCatch("useSourceEncode", 0);
        gradientDomain = config.getAndCatch("gradientDomain", 0);
        useRandomSource = config.getAndCatch("useRandomSource", 0);  
        gradientKernel = config.getAndCatch("gradientKernel", 0); 
        decomposition = config.getAndCatch("decomposition", 0); 
        snapType = config.getAndCatch("snapType", 0);
        breakLoopType = config.get<IndexType>("breakLoopType");
        exchangeStrategy = config.get<IndexType>("exchangeStrategy");
        
        derivatives = ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension);
        solver = ForwardSolver::Factory<ValueType>::Create(dimension, equationType);
        wavefields = Wavefields::Factory<ValueType>::Create(dimension, equationType);
        wavefieldsTemp = Wavefields::Factory<ValueType>::Create(dimension, equationType);
        wavefieldsInversion = Wavefields::Factory<ValueType>::Create(dimension, equationType);
        modelPriori = Modelparameter::Factory<ValueType>::Create(equationType);
        modelPerShot = Modelparameter::Factory<ValueType>::Create(equationType);
        gradient = Gradient::Factory<ValueType>::Create(equationType);
        gradientPerShot = Gradient::Factory<ValueType>::Create(equationType);
        stabilizingFunctionalGradient = Gradient::Factory<ValueType>::Create(equationType);
        gradientOptimization = Optimization::Factory<ValueType>::Create(optimizationType);
    }
}

/*! \brief Print the Configuration
 \param commAll CommunicatorPtr
 \param config Configuration
 \param inversionType inversionType
 \param equationInd equationInd
 */
template <typename ValueType>
void KITGPI::InversionSingle<ValueType>::printConfig(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config, IndexType inversionType, IndexType equationInd)
{
    if (inversionType != 0) {
        HOST_PRINT(commAll, "\n WAVE-Inversion " << dimension << " " << equationType << " " << equationInd << " - LAMA Version\n");
        HOST_PRINT(commAll, "","  - Running on " << commAll->getSize() << " mpi processes -\n\n");    
        if (commAll->getRank() == MASTERGPI) {
            config.print();
        }
        receivers.getModelPerShotSize(commAll, config, NXPerShot, numShotPerSuperShot);
    }
}

/*! \brief Initialization of the InversionSingle
 \param commAll CommunicatorPtr
 \param config Configuration
 \param inversionType inversionType
 \param equationInd equationInd
 \param modelCoordinates modelCoordinates
 \param modelCoordinatesBig modelCoordinatesBig
 \param ctx context
 \param dist dist
 \param distBig distBig
 \param derivativesInversion derivativesInversion
 \param maxiterations maxiterations
 \param model model
 \param workflow workflow
 \param crossGradientDerivative crossGradientDerivative
 \param SLsearch SLsearch
 */
template <typename ValueType>
void KITGPI::InversionSingle<ValueType>::init(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config, IndexType inversionType, IndexType equationInd, Acquisition::Coordinates<ValueType> &modelCoordinates, Acquisition::Coordinates<ValueType> &modelCoordinatesBig, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr &dist, scai::dmemo::DistributionPtr &distBig, typename ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr &derivativesInversion, IndexType maxiterations, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &model, KITGPI::Workflow::Workflow<ValueType> &workflow, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivative, KITGPI::StepLengthSearch<ValueType> &SLsearch)
{
    if (inversionType != 0) {
        /* --------------------------------------- */
        /* coordinate mapping (3D<->1D)            */
        /* --------------------------------------- */
        modelCoordinates.init(config, 1, NXPerShot);
        if (config.get<bool>("useVariableGrid")) {
            CheckParameter::checkVariableGrid(config, commAll, modelCoordinates);
            for (int layer=0;layer<modelCoordinates.getNumLayers();layer++){
                HOST_PRINT(commAll, "\n Number of gridpoints in layer: " << layer << " = " << modelCoordinates.getNGridpoints(layer));
            }
            auto numGridpointsRegular=config.get<IndexType>("NX")*config.get<IndexType>("NY")*config.get<IndexType>("NZ");
            HOST_PRINT(commAll, "\n Number of gripoints total: " << modelCoordinates.getNGridpoints());
            HOST_PRINT(commAll, "\n Percentage of gridpoints of the underlying regular grid given by NX*NY*NZ: "  << (float) modelCoordinates.getNGridpoints()/numGridpointsRegular * 100 << "% \n\n");
        }
        
        /* --------------------------------------- */
        /*          CommunicatorPtr                */
        /* --------------------------------------- */  
        shotDomain = KITGPI::Partitioning::getShotDomain(config, commAll); 

        // Build subsets of processors for the shots
        dmemo::CommunicatorPtr commShot = commAll->split(shotDomain);
        dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());
        SCAI_DMEMO_TASK(commShot)
        
        /* --------------------------------------- */
        /* Distribution                            */
        /* --------------------------------------- */ 
        if (config.get<IndexType>("partitioning") == 0 || config.get<IndexType>("partitioning") == 2) {
            //Block distribution = starting distribution for graph partitioner
            dist = std::make_shared<dmemo::BlockDistribution>(modelCoordinates.getNGridpoints(), commShot);
        } else if (config.get<IndexType>("partitioning") == 1) {
            SCAI_ASSERT(!config.get<bool>("useVariableGrid"), "Grid distribution is not available for the variable grid");
            dist = KITGPI::Partitioning::gridPartition<ValueType>(config, commShot, NXPerShot);
            distInversion = KITGPI::Partitioning::gridPartitionInversion<ValueType>(config, commShot, NXPerShot);
        } else {
            COMMON_THROWEXCEPTION("unknown partitioning method");
        }
        if (config.get<bool>("writeCoordinate") && shotDomain == 0) {
            // every shotdomain owns the same coordinates
            modelCoordinates.writeCoordinates(dist, ctx, config.get<std::string>("coordinateFilename"), config.get<IndexType>("FileFormat"));
        }
        if (useStreamConfig) {
            SCAI_ASSERT_ERROR(config.get<IndexType>("partitioning") == 1, "partitioning != 1 when useStreamConfig"); 
            modelCoordinatesBig.init(configBig);
            distBig = KITGPI::Partitioning::gridPartition<ValueType>(configBig, commShot);
        }  
        numShotDomains = config.get<IndexType>("NumShotDomains"); // total number of shot domains
        Common::checkNumShotDomains(numShotDomains, commAll);
        
        /* --------------------------------------- */
        /* Call partitioner                        */
        /* --------------------------------------- */
        if (config.get<IndexType>("partitioning") == 2) {
            start_t = common::Walltime::get();
            dist = KITGPI::Partitioning::graphPartition(config, ctx, commShot, dist, *derivatives, modelCoordinates);
            end_t = common::Walltime::get();
            HOST_PRINT(commAll, "", "Finished graph partitioning in " << end_t - start_t << " sec.\n\n");
        }    
        
        /* --------------------------------------- */
        /* Calculate derivative matrices           */
        /* --------------------------------------- */
        start_t = common::Walltime::get();
        derivatives->init(dist, ctx, config, modelCoordinates, commShot);
        if (!useStreamConfig) {
            derivativesInversion->init(dist, ctx, config, modelCoordinates, commShot);
        } else {
            derivativesInversion->init(distBig, ctx, config, modelCoordinatesBig, commShot);
        }
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "", "Finished initializing matrices in " << end_t - start_t << " sec.\n\n");
        
        /* --------------------------------------- */
        /* Acquisition geometry                    */
        /* --------------------------------------- */
        ValueType shotIncr = config.getAndCatch("shotIncr", 0.0);
        sources.getAcquisitionSettings(config, shotIncr);
        if (!useStreamConfig) {
            sourceSettings = sources.getSourceSettings(); 
        } else {
            std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;
            sourceSettingsBig = sources.getSourceSettings(); 
            Acquisition::getCutCoord(config, cutCoordinates, sourceSettingsBig, modelCoordinates, modelCoordinatesBig);
            Acquisition::getSettingsPerShot<ValueType>(sourceSettings, sourceSettingsBig, cutCoordinates, modelCoordinates, config.get<IndexType>("BoundaryWidth"));
            sources.setSourceSettings(sourceSettings); // for StepLengthSearch and useSourceEncode
        }
        CheckParameter::checkSources(sourceSettings, modelCoordinates, commShot);
    
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
        if (useSourceEncode == 0) {
            numshots = uniqueShotNos.size();
        } else {
            numshots = numShotDomains;
        }
        SCAI_ASSERT_ERROR(numshots >= numShotDomains, "numshots = " + std::to_string(numshots) + ", numShotDomains = " + std::to_string(numShotDomains));
        if (useSourceEncode == 0 && useRandomSource == 0) {
            SCAI_ASSERT_ERROR(numshots % numShotDomains == 0, "numshots = " + std::to_string(numshots) + ", numShotDomains = " + std::to_string(numShotDomains));
        }    
        sources.writeShotIndsIncr(commAll, config, uniqueShotNos);
        Acquisition::writeCutCoordToFile(commAll, config, cutCoordinates, uniqueShotNos, NXPerShot);
        
        if (useRandomSource != 0) {  
            shotDist = dmemo::blockDistribution(numShotDomains, commInterShot);
            maxcount = ceil((ValueType)maxiterations * numShotDomains / numshots * 1.2);
        } else {
            shotDist = dmemo::blockDistribution(numshots, commInterShot);
        }
        if (config.get<IndexType>("useReceiversPerShot") == 0) {
            receivers.init(config, modelCoordinates, ctx, dist);
        }
        
        if (uniqueShotNos.size() == sourceSettings.size() && uniqueShotNos.size() > 1) {
            if (config.get<bool>("writeInvertedSource"))
                sources.getSeismogramHandler().allocateCOP(numshots, tStepEnd);
            if (receivers.getNumTracesGlobal() == numShotPerSuperShot)
                receivers.getSeismogramHandler().allocateCOP(numshots, tStepEnd);
        }
        
        /* --------------------------------------- */
        /* Wavefields                              */
        /* --------------------------------------- */ 
        wavefields->init(ctx, dist, numRelaxationMechanisms);
        wavefieldsInversion->init(ctx, distInversion, numRelaxationMechanisms);
        if ((gradientKernel == 2 || gradientKernel == 3) && decomposition == 0)
            wavefieldsTemp->init(ctx, dist, numRelaxationMechanisms);
        if ((gradientKernel == 2 || gradientKernel == 3) && decomposition != 0)
            snapType = decomposition + 3;
        
        /* --------------------------------------- */
        /* Modelparameter                          */
        /* --------------------------------------- */
        model->prepareForInversion(config, commShot);
        modelPriori->prepareForInversion(config, commShot);
        modelPerShot->prepareForInversion(config, commShot); // prepareForInversion is necessary for modelPerShot to calculate gradient.
        if (!useStreamConfig) {    
            model->init(config, ctx, dist, modelCoordinates);
            modelPriori->init(config, ctx, dist, modelCoordinates);
        } else { 
            model->init(configBig, ctx, distBig, modelCoordinatesBig);   
            modelPriori->init(configBig, ctx, distBig, modelCoordinatesBig);
        }
        
        /* --------------------------------------- */
        /* Forward solver                          */
        /* --------------------------------------- */
        if (!useStreamConfig) {
            solver->initForwardSolver(config, *derivatives, *wavefields, *model, modelCoordinates, ctx, config.get<ValueType>("DT"));
        }
        
        /* --------------------------------------- */
        /* True data and Adjoint sources           */
        /* --------------------------------------- */
        if (config.get<IndexType>("useReceiversPerShot") == 0) {
            receiversTrue.init(config, modelCoordinates, ctx, dist);
            receiversStart.init(config, modelCoordinates, ctx, dist);
            adjointSources.init(config, modelCoordinates, ctx, dist);
        }
        if (receivers.getNumTracesGlobal() == numShotPerSuperShot) {
            receiversTrue.getSeismogramHandler().allocateCOP(numshots, tStepEnd);
            receiversStart.getSeismogramHandler().allocateCOP(numshots, tStepEnd);
            adjointSources.getSeismogramHandler().allocateCOP(numshots, tStepEnd);
        }
        seismogramTaper1D.init(std::make_shared<dmemo::NoDistribution>(tStepEnd), ctx, 1);
        
        /* --------------------------------------- */
        /* Misfit                                  */
        /* --------------------------------------- */
        misfitPerIt = scai::lama::fill<scai::lama::DenseVector<ValueType>>(numshots, 0);
        
        /* --------------------------------------- */
        /* Step length search                      */
        /* --------------------------------------- */
        SLsearch.initLogFile(commAll, logFilename, misfitType, config.getAndCatch("steplengthType", 2), workflow.getInvertForParameters().size(), config.getAndCatch("saveCrossGradientMisfit", 0));
        
        /* --------------------------------------- */
        /* Source estimation                       */
        /* --------------------------------------- */
        // calculate source dist
        lama::DenseVector<IndexType> sourcecoords = getsourcecoordinates(sourceSettings, modelCoordinates);
        dmemo::DistributionPtr dist_sources = Acquisition::calcDistribution(sourcecoords, dist);
        if (config.get<IndexType>("useSourceSignalInversion") != 0)
            sourceEst.init(config, workflow, ctx, dist_sources, sourceSignalTaper);
        
        /* --------------------------------------- */
        /* Frequency filter                        */
        /* --------------------------------------- */
        if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0)
            freqFilter.init(config.get<ValueType>("DT"), tStepEnd);
        
        /* --------------------------------------- */
        /* Gradients                               */
        /* --------------------------------------- */
        if (!useStreamConfig) {
            gradient->init(ctx, dist);      
            stabilizingFunctionalGradient->init(ctx, dist);
            crossGradientDerivative->init(ctx, dist);  
        } else {
            gradient->init(ctx, distBig);
            stabilizingFunctionalGradient->init(ctx, distBig);
            crossGradientDerivative->init(ctx, distBig);
        }
        gradientPerShot->init(ctx, dist); 
        gradient->calcGaussianKernel(commAll, *model, config);
        
        gradient->prepareForInversion(config);
        gradientPerShot->prepareForInversion(config);  
        stabilizingFunctionalGradient->prepareForInversion(config);
        crossGradientDerivative->prepareForInversion(config);  
        
        /* --------------------------------------- */
        /* Gradient taper                          */
        /* --------------------------------------- */
        if (config.get<bool>("useGradientTaper")) {
            if (!useStreamConfig) {
                gradientTaper1D.init(dist, ctx, 1);
                gradientTaper1D.read(config.get<std::string>("gradientTaperName"), config.get<IndexType>("FileFormat"));
            } else {
                gradientTaper1D.init(distBig, ctx, 1);
                gradientTaper1D.read(configBig.get<std::string>("gradientTaperName"), config.get<IndexType>("FileFormat"));
            }  
        }
        wavefieldTaper2D.initAverageMatrix(config, distInversion, dist, ctx); 
        Acquisition::Coordinates<ValueType> modelCoordinatesInversion(config, config.getAndCatch("DHInversion", 1), NXPerShot);
        wavefieldTaper2D.calcAverageMatrix(modelCoordinates, modelCoordinatesInversion);
        
        /* --------------------------------------- */
        /* Gradient preconditioning                */
        /* --------------------------------------- */
        energyPrecond.init(distInversion, config);
        energyPrecondReflect.init(distInversion, config);
        
        /* --------------------------------------- */
        /* Gradient optimization                   */
        /* --------------------------------------- */
        gradientOptimization->init(dist);
        
        /* --------------------------------------- */
        /*       Modelparameter preparation        */
        /* --------------------------------------- */    
        if (inversionType == 3 || model->getParameterisation() == 2 || model->getParameterisation() == 1) {
            // in case of using porosity and saturation in inversion
            model->calcRockMatrixParameter(config);     
        }
    }
}

/*! \brief Estimate memory of the InversionSingle
 \param commAll CommunicatorPtr
 \param config Configuration
 \param inversionType inversionType
 \param equationInd equationInd
 \param dist dist
 \param modelCoordinates modelCoordinates
 \param model model
 */
template <typename ValueType>
void KITGPI::InversionSingle<ValueType>::estimateMemory(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config, IndexType inversionType, IndexType equationInd, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> modelCoordinates, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model)
{
    if (inversionType != 0) {
        ValueType memDerivatives = derivatives->estimateMemory(config, dist, modelCoordinates);
        ValueType memWavefileds = wavefields->estimateMemory(dist, numRelaxationMechanisms);
        IndexType NT = tStepEnd;
        if (useSourceEncode == 0 && (gradientDomain == 1 || gradientDomain == 2)) {            
            NT *= 2;
        } else if (gradientDomain == 3) {            
            NT = numShotPerSuperShot * 2;
        }
        if (dimension.compare("3d") == 0) {
            memWavefiledsStorage = memWavefileds * NT / pow(config.getAndCatch("DHInversion", 1), 3);
        } else {
            memWavefiledsStorage = memWavefileds * NT / pow(config.getAndCatch("DHInversion", 1), 2);
        }
        ValueType memModel = model->estimateMemory(dist);
        ValueType memSolver = solver->estimateMemory(config, dist, modelCoordinates);
        ValueType memTotal = memDerivatives + memWavefileds + memModel + memSolver + memWavefiledsStorage;

        HOST_PRINT(commAll, "============== Memory Estimation " << equationType << " " << equationInd << ": ===============\n\n")
        HOST_PRINT(commAll, " -  Derivative Matrices \t" << memDerivatives << " MB\n");
        HOST_PRINT(commAll, " -  Wavefield vectors \t\t" << memWavefileds << " MB\n");
        HOST_PRINT(commAll, " -  Forward wavefield storage \t" << memWavefiledsStorage << " MB\n");
        HOST_PRINT(commAll, " -  Model Vectors \t\t" << memModel << " MB\n");
        HOST_PRINT(commAll, " -  Boundary Condition Vectors \t" << memSolver << " MB\n");
        HOST_PRINT(commAll, "\n Memory Usage (total / per partition): \n " << memTotal << " / " << memTotal / dist->getNumPartitions() << " MB ");
        if (numShotDomains > 1)
            HOST_PRINT(commAll, "\n Total Memory Usage (" << numShotDomains << " shot domains): \n " << memTotal * numShotDomains << " MB  ");

        HOST_PRINT(commAll, "\n\n===================================================\n")  
    }
}

/*! \brief Initialization of the InversionSingle stage
 \param commAll CommunicatorPtr
 \param config Configuration
 \param inversionType inversionType
 \param equationInd equationInd
 \param ctx context
 \param workflow workflow
 \param dataMisfit dataMisfit
 \param breakLoop breakLoop
 \param dist dist
 \param breakLoopEM breakLoopEM
 */
template <typename ValueType>
void KITGPI::InversionSingle<ValueType>::initStage(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config, IndexType inversionType, IndexType equationInd, scai::hmemo::ContextPtr ctx, KITGPI::Workflow::Workflow<ValueType> &workflow, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr dataMisfit, bool breakLoop, scai::dmemo::DistributionPtr dist, bool &breakLoopEM) 
{
    if (inversionType != 0 && (breakLoop == false || breakLoopType == 2)) {
        if (equationInd == 1 && (exchangeStrategy == 4 || exchangeStrategy == 6))
            breakLoopEM = true;

        HOST_PRINT(commAll, "\nChange workflow stage " << equationType << " " << equationInd << "\n");
        workflow.changeStage(config, *dataMisfit, steplengthInit);
                
        workflow.printParameters(commAll);

        gradientCalculation.allocate(config, dist, distInversion, ctx, workflow, numShotPerSuperShot);
        seismogramTaper1D.calcTimeDampingTaper(workflow.getTimeDampingFactor(), config.get<ValueType>("DT"));  

        if (workflow.getLowerCornerFreq() != 0.0 && workflow.getUpperCornerFreq() != 0.0)
            freqFilter.calc(transFcnFmly, "bp", workflow.getFilterOrder(), workflow.getLowerCornerFreq(), workflow.getUpperCornerFreq());
        else if (workflow.getLowerCornerFreq() != 0.0 && workflow.getUpperCornerFreq() == 0.0)
            freqFilter.calc(transFcnFmly, "lp", workflow.getFilterOrder(), workflow.getLowerCornerFreq());
        else if (workflow.getLowerCornerFreq() == 0.0 && workflow.getUpperCornerFreq() != 0.0)
            freqFilter.calc(transFcnFmly, "hp", workflow.getFilterOrder(), workflow.getUpperCornerFreq()); 
        
        if (workflow.skipDT > 1) {
            HOST_PRINT(commAll, "\nForward wavefield storage reduces from " << memWavefiledsStorage << " MB to " << memWavefiledsStorage / workflow.skipDT << " MB (skipDT = " << workflow.skipDT << ")\n");
        }
        IndexType NT = tStepEnd;
        if (gradientDomain != 0) {
            NT = 1;
        }
        wavefieldrecord.clear();
        for (IndexType tStep = 0; tStep < NT; tStep++) {
            if (tStep % workflow.skipDT == 0) {
                // Only the local shared_ptr can be used to initialize a std::vector
                wavefieldPtr wavefieldsInversion = Wavefields::Factory<ValueType>::Create(dimension, equationType);
                wavefieldsInversion->init(ctx, distInversion, numRelaxationMechanisms);
                wavefieldrecord.push_back(wavefieldsInversion);
            }
        }
        if ((gradientKernel == 2 || gradientKernel == 3) && decomposition == 0) {
            wavefieldrecordReflect.clear();
            for (IndexType tStep = 0; tStep < NT; tStep++) {
                if (tStep % workflow.skipDT == 0) {
                    wavefieldPtr wavefieldsInversion = Wavefields::Factory<ValueType>::Create(dimension, equationType);
                    wavefieldsInversion->init(ctx, distInversion, numRelaxationMechanisms);
                    wavefieldrecordReflect.push_back(wavefieldsInversion);
                }
            }
        }
        std::vector<IndexType> temp(numshots, 0);
        shotHistory = temp;
        std::vector<IndexType> temp2(misfitType.length() - 2, 0);
        misfitTypeHistory = temp2;
    }
}
    
/*! \brief Calculate gradient of the InversionSingle
 \param commAll CommunicatorPtr
 \param dist dist
 \param model model
 \param config config
 \param modelCoordinates modelCoordinates
 \param modelCoordinatesBig modelCoordinatesBig
 \param workflow workflow
 \param dataMisfit dataMisfit
 \param crossGradientDerivative crossGradientDerivative
 \param SLsearch SLsearch
 \param modelTaper2DJoint modelTaper2DJoint
 \param maxiterations maxiterations
 \param useRTM useRTM
 \param breakLoop breakLoop
 \param ctx context
 \param seedtime seedtime
 \param inversionType inversionType
 \param equationInd equationInd
 \param breakLoopEM breakLoopEM
 \param modelEM modelEM
 \param configEM configEM
 \param modelCoordinatesEM modelCoordinatesEM
 \param workflowEM workflowEM
 \param dataMisfitEM dataMisfitEM
 \param crossGradientDerivativeEM crossGradientDerivativeEM
 \param derivativesInversionEM derivativesInversionEM
 */
template <typename ValueType>
void KITGPI::InversionSingle<ValueType>::calcGradient(scai::dmemo::CommunicatorPtr commAll, scai::dmemo::DistributionPtr dist, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &model, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, KITGPI::Workflow::Workflow<ValueType> &workflow, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfit, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivative, KITGPI::StepLengthSearch<ValueType> &SLsearch, Taper::Taper2D<ValueType> modelTaper2DJoint, IndexType maxiterations, IndexType &useRTM, bool &breakLoop, scai::hmemo::ContextPtr ctx, IndexType &seedtime, IndexType inversionType, IndexType equationInd, bool &breakLoopEM, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &modelEM, KITGPI::Configuration::Configuration configEM, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesEM, KITGPI::Workflow::Workflow<ValueType> &workflowEM, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfitEM, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivativeEM, typename ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr &derivativesInversionEM)
{          
    if (inversionType != 0 && (breakLoop == false || breakLoopType == 2 || useRTM != 0)) {
        IndexType shotNumber;
        IndexType shotIndTrue = 0; 
        IndexType shotIndIncr = 0;  
        
        dmemo::CommunicatorPtr commShot = dist->getCommunicatorPtr();
        dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());
        SCAI_DMEMO_TASK(commShot)
        
        if (useRTM != 0)
            HOST_PRINT(commAll, "\nStart inverse time migration after stage " << workflow.workflowStage + 1 << " in "<< equationType << " " << equationInd << "");
        if (equationInd == 1 && (exchangeStrategy == 3 || exchangeStrategy == 5))
            breakLoopEM = true;
        
        // Begin of one seismic model update 
        HOST_PRINT(commAll, "\n=================================================");
        HOST_PRINT(commAll, "\n=========== Workflow stage " << workflow.workflowStage + 1 << " of " << workflow.maxStage << " ===============");
        HOST_PRINT(commAll, "\n============     Iteration " << workflow.iteration + 1 << "       ==============");
        HOST_PRINT(commAll, "\n=================== " << equationType << " " << equationInd << " =======================\n\n");
        start_t = common::Walltime::get();
        
        if (!useStreamConfig) { 
            /* Update model for fd simulation (averaging, inverse Density ...) */
            *modelPerShot = *model;
            modelPerShot->prepareForModelling(modelCoordinates, ctx, dist, commShot); 
            solver->prepareForModelling(*modelPerShot, config.get<ValueType>("DT"));
        }
        
        sources.calcUniqueShotInds(commAll, config, shotHistory, maxcount, seedtime);
        std::vector<IndexType> uniqueShotInds = sources.getUniqueShotInds();
        std::vector<IndexType> shotIndsIncr = sources.getShotIndsIncr();
        Acquisition::writeRandomShotNosToFile(commAll, logFilename, uniqueShotNos, uniqueShotInds, workflow.workflowStage + 1, workflow.iteration, useRandomSource);
        dataMisfit->init(config, misfitTypeHistory, numshots, useRTM, model->getVmin(), seedtime); // in case of that random misfit function is used
        sources.calcSourceSettingsEncode(commAll, config, seedtime, workflow.getLowerCornerFreq(), workflow.getUpperCornerFreq()); // for sourceFC
        if (useSourceEncode != 0) {
            sourceSettingsEncode = sources.getSourceSettingsEncode();
            Acquisition::calcuniqueShotNo(uniqueShotNosEncode, sourceSettingsEncode);
        }
        sources.writeSourceFC(commAll, config, workflow.workflowStage + 1, workflow.iteration);
        sources.writeSourceEncode(commAll, config, workflow.workflowStage + 1, workflow.iteration);
            
        if (workflow.iteration == 0 && commInterShot->getRank() == 0) {
            /* only shot domain 0 writes output */
            model->write((config.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration)), config.get<IndexType>("FileFormat"));
        }
        
        /* --------------------------------------- */
        /*        Loop over Seismic shots          */
        /* --------------------------------------- */
        gradient->resetGradient(); // reset gradient because gradient is a sum of all gradientsPerShot gradients+=gradientPerShot
        crossGradientDerivative->resetGradient();
        misfitPerIt = 0;
        
        IndexType gradientKernelPerIt = 0; 
        if (gradientKernel == 3) {
            IndexType numSwitch = gradientKernel - 2;  
            if ((workflow.iteration / numSwitch) % 2 == 0) {
                gradientKernelPerIt = 1;
            } else {
                gradientKernelPerIt = 2;
            }            
            if (decomposition == 0 && workflow.iteration % (numSwitch * 2) == 0) {
                model->resetReflectivity();
            }
        } else {
            gradientKernelPerIt = gradientKernel;
        }
    
        std::vector<bool> invertForParameters = workflow.getInvertForParameters();
        if (gradientKernelPerIt == 1 && decomposition == 0) { // the last element of invertForParameters is invertForReflectivity
            invertForParameters[invertForParameters.size()-1] = true;
//                     if (!useStreamConfig) {
//                         model->calcReflectivity(modelCoordinates, *derivatives, config.get<ValueType>("DT"));
//                     } else {
//                         modelPerShot->calcReflectivity(modelCoordinates, *derivatives, config.get<ValueType>("DT"));
//                     }
        }
        gradient->setInvertForParameters(invertForParameters);
        gradientPerShot->setInvertForParameters(invertForParameters);
        crossGradientDerivative->setInvertForParameters(invertForParameters);
        stabilizingFunctionalGradient->setInvertForParameters(invertForParameters);
        workflow.setInvertForParameters(invertForParameters);
        
        IndexType localShotInd = 0;     
        for (IndexType shotInd = shotDist->lb(); shotInd < shotDist->ub(); shotInd++) {
            shotIndTrue = uniqueShotInds[shotInd];
            shotIndIncr = shotIndsIncr[shotInd]; // it is not compatible with useSourceEncode != 0
            localShotInd++;
            
            std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
            if (useSourceEncode == 0) {
                shotNumber = uniqueShotNos[shotIndTrue];
                Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
            } else {
                shotNumber = uniqueShotNosEncode[shotIndTrue];
                Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettingsEncode, shotNumber);
            }                    
            sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);
            if (uniqueShotNos.size() == sourceSettings.size() && uniqueShotNos.size() > 1) {
                sources.getSeismogramHandler().setShotInd(shotIndTrue, shotIndIncr);
            }

            IndexType shotIndPerShot = shotIndTrue;
            if (useSourceEncode == 3) {
                Acquisition::getuniqueShotInd(shotIndPerShot, sourceSettingsEncode, shotNumber);
            }
            if (useStreamConfig) {
                HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Switch to model subset\n");
                
                model->getModelPerShot(*modelPerShot, dist, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotIndPerShot)); 
                modelPerShot->prepareForModelling(modelCoordinates, ctx, dist, commShot); 
                solver->initForwardSolver(config, *derivatives, *wavefields, *modelPerShot, modelCoordinates, ctx, config.get<ValueType>("DT"));
                solver->prepareForModelling(*modelPerShot, config.get<ValueType>("DT"));
            }
            CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(config, sourceSettingsShot, *modelPerShot, modelCoordinates, shotNumber);
            HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "), local shot " << localShotInd << " of " << shotDist->getLocalSize() << ": Started\n");
                
            if (config.get<IndexType>("useReceiversPerShot") != 0) {
                receivers.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncode);
                receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncode);
                adjointSources.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncode);
            }
            if (uniqueShotNos.size() == sourceSettings.size() && uniqueShotNos.size() > 1 && receivers.getNumTracesGlobal() == numShotPerSuperShot) {
                receivers.getSeismogramHandler().setShotInd(shotIndTrue, shotIndIncr);
                receiversTrue.getSeismogramHandler().setShotInd(shotIndTrue, shotIndIncr);
                adjointSources.getSeismogramHandler().setShotInd(shotIndTrue, shotIndIncr);
            }
            /* Read field data (or pseudo-observed data, respectively) */
            if (useSourceEncode == 0) {
                receiversTrue.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".shot_" + std::to_string(shotNumber), 1);
            } else {
                receiversTrue.encode(config, config.get<std::string>("fieldSeisName"), shotNumber, sourceSettingsEncode, 1);
            }
                                
            if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0){
                receiversTrue.getSeismogramHandler().filter(freqFilter);
            }
            
            if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0)
                sources.getSeismogramHandler().filter(freqFilter);
            
            /* Source time function inversion */                    
            std::string filenameSyn = config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration);
            std::string filenameObs = config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration);
            
            if (workflow.getMinOffset() != 0 || workflow.getMaxOffset() != 0 || misfitType.compare("l3") == 0 || multiMisfitType.find('3') != std::string::npos || misfitType.compare("l4") == 0 || multiMisfitType.find('4') != std::string::npos){
                if (useSourceEncode == 0) {
                    sourceEst.calcOffsetMutes(sources, receiversTrue, workflow.getMinOffset(), workflow.getMaxOffset(), shotIndTrue, modelCoordinates);
                    sourceEst.applyOffsetMute(config, shotIndTrue, receiversTrue);
                } else {
                    sourceEst.calcOffsetMutesEncode(commShot, shotNumber, config, workflow, modelCoordinates, ctx, dist, sourceSettingsEncode, receiversTrue);
                    receiversTrue.decode(config, filenameObs, shotNumber, sourceSettingsEncode, 0);
                    sourceEst.applyOffsetMuteEncode(commShot, shotNumber, config, sourceSettingsEncode, receiversTrue);
                }
                if (config.get<IndexType>("useSourceSignalInversion") == 2 || misfitType.compare("l3") == 0 || multiMisfitType.find('3') != std::string::npos) {
                    if (config.get<IndexType>("useSourceSignalTaper") == 2) {
                        sourceSignalTaper.calcCosineTaper(sources.getSeismogramHandler(), workflow.getLowerCornerFreq(), workflow.getUpperCornerFreq(), config.get<ValueType>("DT"), ctx);
                    }
                    if (useSourceEncode == 0) {
                        sourceEst.calcRefTraces(config, shotIndTrue, receiversTrue, sourceSignalTaper);
                    } else {
                        sourceEst.calcRefTracesEncode(commShot, shotNumber, config, modelCoordinates, ctx, dist, sourceSettingsEncode, receiversTrue, sourceSignalTaper);
                    }
                    sourceEst.setRefTracesToSource(sources, receiversTrue, sourceSettingsEncode, shotIndTrue, shotNumber);
                }
            }
            if (config.get<IndexType>("useSourceSignalInversion") != 0){
                if (workflow.iteration == 0 || shotHistory[shotIndTrue] != 0 || useSourceEncode != 0) {
                    HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "), local shot " << localShotInd << " of " << shotDist->getLocalSize() << ": Source Time Function Inversion\n");

                    wavefields->resetWavefields();

                    for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                        solver->run(receivers, sources, *modelPerShot, *wavefields, *derivatives, tStep);
                    }
                    solver->resetCPML();
                    
                    /* Normalize observed and synthetic data */
                    if (config.get<IndexType>("normalizeTraces") == 3 || misfitType.compare("l6") == 0 || multiMisfitType.find('6') != std::string::npos) {
                        ValueType frequencyAGC = config.get<ValueType>("CenterFrequencyCPML");
                        if (workflow.getUpperCornerFreq() != 0.0) {
                            frequencyAGC = (workflow.getLowerCornerFreq() + workflow.getUpperCornerFreq()) / 2;
                        }
                        receiversTrue.getSeismogramHandler().setFrequencyAGC(frequencyAGC);       
                        receiversTrue.getSeismogramHandler().calcInverseAGC();
                        receiversTrue.getSeismogramHandler().write(5, config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);    
                        receivers.getSeismogramHandler().setInverseAGC(receiversTrue.getSeismogramHandler());
                    }
                    receivers.getSeismogramHandler().normalize(config.get<IndexType>("normalizeTraces"));
                    receiversTrue.getSeismogramHandler().normalize(config.get<IndexType>("normalizeTraces"));
                    
                    if (useSourceEncode == 0) {
                        sourceEst.applyOffsetMute(config, shotIndTrue, receivers);
                        sourceEst.estimateSourceSignal(receivers, receiversTrue, shotIndTrue, shotNumber);
                    } else {
                        receivers.decode(config, filenameSyn, shotNumber, sourceSettingsEncode, 0);
                        receiversTrue.decode(config, filenameObs, shotNumber, sourceSettingsEncode, 0);
                        sourceEst.applyOffsetMuteEncode(commShot, shotNumber, config, sourceSettingsEncode, receivers);
                        sourceEst.estimateSourceSignalEncode(commShot, shotNumber, config, modelCoordinates, ctx, dist, sourceSettingsEncode, receivers, receiversTrue);
                    }
                }
                if (useSourceEncode == 0) {
                    sourceEst.applyFilter(sources, shotNumber, sourceSettings);
                } else {
                    sourceEst.applyFilter(sources, shotNumber, sourceSettingsEncode);
                }
                
                if (config.get<IndexType>("useSourceSignalTaper") != 0) {
                    if (config.get<IndexType>("useSourceSignalTaper") == 2) {
                        sourceSignalTaper.calcCosineTaper(sources.getSeismogramHandler(), workflow.getLowerCornerFreq(), workflow.getUpperCornerFreq(), config.get<ValueType>("DT"), ctx);
                    }
                    sourceSignalTaper.apply(sources.getSeismogramHandler());
                }
            }
            if (config.get<bool>("writeInvertedSource") && (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1))
                sources.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("sourceSeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
            
            if (config.get<IndexType>("useSeismogramTaper") > 1 && config.get<IndexType>("useSeismogramTaper") != 5) {
                seismogramTaper2D.init(receiversTrue.getSeismogramHandler());
                if (config.get<IndexType>("useSeismogramTaper") == 4) {
                    seismogramTaper2D.read(config.get<std::string>("seismogramTaperName") + ".misfitCalc.shot_" + std::to_string(shotNumber) + ".mtx");
                } else {
                    seismogramTaper2D.read(config.get<std::string>("seismogramTaperName") + ".shot_" + std::to_string(shotNumber) + ".mtx");
                }
                seismogramTaper2D.apply(receiversTrue.getSeismogramHandler()); 
            }
            if (config.get<IndexType>("useSeismogramTaper") == 5 && (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1)) {
                seismogramTaper1D.calcCosineTaper(receiversTrue.getSeismogramHandler(), workflow.getUpperCornerFreq(), workflow.getUpperCornerFreq(), config.get<ValueType>("DT"), ctx, 1);
            } 
            seismogramTaper1D.apply(receiversTrue.getSeismogramHandler());                                        
            
            /* Normalize observed and synthetic data */
            if (config.get<IndexType>("normalizeTraces") == 3 || misfitType.compare("l6") == 0 || multiMisfitType.find('6') != std::string::npos) {
                if (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1) {
                    if (config.get<IndexType>("useSourceSignalInversion") == 0) {
                        ValueType frequencyAGC = config.get<ValueType>("CenterFrequencyCPML");
                        if (workflow.getUpperCornerFreq() != 0.0) {
                            frequencyAGC = (workflow.getLowerCornerFreq() + workflow.getUpperCornerFreq()) / 2;
                        }
                        receiversTrue.getSeismogramHandler().setFrequencyAGC(frequencyAGC);       
                        receiversTrue.getSeismogramHandler().calcInverseAGC();
                        // to write inverseAGC matrix.
                        receiversTrue.getSeismogramHandler().write(5, config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);          
                    }
                } else {
                    // to read inverseAGC matrix.
                    receiversTrue.getSeismogramHandler().read(5, config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber));             
                }
                receivers.getSeismogramHandler().setInverseAGC(receiversTrue.getSeismogramHandler());
            }
            receiversTrue.getSeismogramHandler().normalize(config.get<IndexType>("normalizeTraces"));
            
            if (useSourceEncode == 0 && (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1)) {
                receiversTrue.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates); 
            } else if (useSourceEncode != 0) {
                if (config.getAndCatch("writeAdjointSource", false)) {
                    receiversTrue.decode(config, filenameObs, shotNumber, sourceSettingsEncode, 1);
                } else {
                    receiversTrue.decode(config, filenameObs, shotNumber, sourceSettingsEncode, 0);
                }
                receiversTrue.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), filenameObs + ".shot_" + std::to_string(shotNumber), modelCoordinates); 
            }
            
            if (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1) {                        
                if (config.get<IndexType>("useReceiversPerShot") != 0) {
                    receiversStart.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncode);
                }
                if (uniqueShotNos.size() == sourceSettings.size() && uniqueShotNos.size() > 1 && receivers.getNumTracesGlobal() == numShotPerSuperShot) {
                    receiversStart.getSeismogramHandler().setShotInd(shotIndTrue, shotIndIncr);
                }
            
                if (useSourceEncode == 0) {
                    receiversStart.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber), 1);
                } else {
                    receiversStart.encode(config, config.get<std::string>("SeismogramFilename"), shotNumber, sourceSettingsEncode, 1);
                }
                if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0){
                    receiversStart.getSeismogramHandler().filter(freqFilter);
                }
                if (config.get<IndexType>("useSeismogramTaper") > 1 && config.get<IndexType>("useSeismogramTaper") != 5) {
                    seismogramTaper2D.apply(receiversStart.getSeismogramHandler()); 
                }
                seismogramTaper1D.apply(receiversStart.getSeismogramHandler());
                
                if (config.get<IndexType>("normalizeTraces") == 3)
                    receiversStart.getSeismogramHandler().setInverseAGC(receiversTrue.getSeismogramHandler());
                
                receiversStart.getSeismogramHandler().normalize(config.get<IndexType>("normalizeTraces"));
                
                if (config.getAndCatch("writeAdjointSource", false))
                    receiversStart.decode(config, config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1), shotNumber, sourceSettingsEncode, 1);
                
                if (useSourceEncode == 0) {
                    sourceEst.applyOffsetMute(config, shotIndTrue, receiversStart);
                } else {
                    sourceEst.applyOffsetMuteEncode(commShot, shotNumber, config, sourceSettingsEncode, receiversStart);
                }
                
                receiversStart.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
            }
            
            /* --------------------------------------- */
            /*        Forward modelling 1              */
            /* --------------------------------------- */
            HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Start time stepping with " << tStepEnd << " time steps\n");

            ValueType DTinv = 1.0 / config.get<ValueType>("DT");
            lama::DenseVector<ValueType> compensation;
            if (gradientKernelPerIt == 2 && decomposition == 0) { 
                HOST_PRINT(commAll, "================ initWholeSpace receivers ===============\n");
                sourcesReflect.initWholeSpace(config, modelCoordinates, ctx, dist, receivers.getSeismogramTypes());
            }   
            typename ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::SourceReceiverImplPtr SourceReceiverReflect(ForwardSolver::SourceReceiverImpl::Factory<ValueType>::Create(dimension, equationType, sources, sourcesReflect, *wavefieldsTemp));

            start_t_shot = common::Walltime::get();
            wavefields->resetWavefields();
            energyPrecond.resetApproxHessian();
        
            for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                *wavefieldsTemp = *wavefields;

                solver->run(receivers, sources, *modelPerShot, *wavefields, *derivatives, tStep);
                
                if ((gradientKernelPerIt == 2 && decomposition == 0) || decomposition != 0) { 
                    //calculate temporal derivative of wavefield
                    *wavefieldsTemp -= *wavefields;
                    *wavefieldsTemp *= -DTinv; // wavefieldsTemp will be gathered by sourcesReflect
                    if (gradientKernelPerIt == 2 && decomposition == 0) 
                        SourceReceiverReflect->gatherSeismogram(tStep);
                    if (decomposition != 0) 
                        wavefields->decompose(decomposition, *wavefieldsTemp, *derivatives);
                }
                if (tStep % workflow.skipDT == 0 && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep >= tStepEnd / 2))) {
                    if (config.getAndCatch("compensation", 0)) {
                        compensation = modelPerShot->getCompensation(config.get<ValueType>("DT"), tStep);
                        *wavefieldsInversion = *wavefields;
                        *wavefieldsInversion *= compensation;
                        if (config.getAndCatch("DHInversion", 1) > 1)
                            wavefieldsInversion->applyTransform(wavefieldTaper2D.getAverageMatrix(), *wavefieldsInversion);
                    } else {
                        if (config.getAndCatch("DHInversion", 1) > 1) {
                            wavefieldsInversion->applyTransform(wavefieldTaper2D.getAverageMatrix(), *wavefields);
                        } else {
                            *wavefieldsInversion = *wavefields;
                        }
                    }
                    if (gradientDomain == 0 || tStep == 0) {
                        *wavefieldrecord[floor(tStep / workflow.skipDT + 0.5)] = *wavefieldsInversion;
                    } 
                    if (gradientDomain != 0 && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep >= tStepEnd / 2))) {
                        gradientCalculation.gatherWavefields(*wavefieldsInversion, sources.getSourceFC(shotIndTrue), workflow, tStep, config.get<ValueType>("DT"));
                    }
                    energyPrecond.intSquaredWavefields(*wavefieldsInversion, config.get<ValueType>("DT"));
                }       
                
                if (workflow.workflowStage == 0 && workflow.iteration == 0 && gradientDomain == 0 && config.getAndCatch("snapType", 0) > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT")) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), config.get<ValueType>("DT")) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT"))) % Common::time2index(config.get<ValueType>("tincSnapshot"), config.get<ValueType>("DT")) == 0) {
                    wavefields->write(snapType, config.getAndCatch<std::string>("WavefieldFileName", "") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) +  + ".shot_" + std::to_string(shotNumber) + ".source", tStep, *derivatives, *modelPerShot, config.get<IndexType>("FileFormat"));
                }
            }
            solver->resetCPML();
            
            if (gradientKernelPerIt == 2 && decomposition == 0) { 
                HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Start reflection forward \n");
                scai::lama::DenseVector<ValueType> reflectivity;
                reflectivity = modelPerShot->getReflectivity();
                dataMisfit->calcReflectSources(sourcesReflect, reflectivity);
                wavefields->resetWavefields(); 
                energyPrecondReflect.resetApproxHessian();
                bool isReflect = true;
                
                for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                    
                    solver->run(adjointSources, sourcesReflect, *modelPerShot, *wavefields, *derivatives, tStep);
                    
                    if (tStep % workflow.skipDT == 0 && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep >= tStepEnd / 2))) {
                        if (config.getAndCatch("compensation", 0)) {
                            compensation = modelPerShot->getCompensation(config.get<ValueType>("DT"), tStep);
                            *wavefieldsInversion = *wavefields;
                            *wavefieldsInversion *= compensation;                                
                            if (config.getAndCatch("DHInversion", 1) > 1)
                                wavefieldsInversion->applyTransform(wavefieldTaper2D.getAverageMatrix(), *wavefieldsInversion);
                        } else {
                            if (config.getAndCatch("DHInversion", 1) > 1) {
                                wavefieldsInversion->applyTransform(wavefieldTaper2D.getAverageMatrix(), *wavefields);
                            } else {
                                *wavefieldsInversion = *wavefields;
                            }
                        }
                        if (gradientDomain == 0 || tStep == 0) {
                            *wavefieldrecordReflect[floor(tStep / workflow.skipDT + 0.5)] = *wavefieldsInversion;
                        } 
                        if (gradientDomain != 0 && (useSourceEncode == 0 || (useSourceEncode != 0 && tStep >= tStepEnd / 2))) {
                            gradientCalculation.gatherWavefields(*wavefieldsInversion, sources.getSourceFC(shotIndTrue), workflow, tStep, config.get<ValueType>("DT"), isReflect);
                        }
                        energyPrecondReflect.intSquaredWavefields(*wavefieldsInversion, config.get<ValueType>("DT"));
                    }
                    if (workflow.workflowStage == 0 && workflow.iteration == 0 && gradientDomain == 0 && config.getAndCatch("snapType", 0) > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT")) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), config.get<ValueType>("DT")) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), config.get<ValueType>("DT"))) % Common::time2index(config.get<ValueType>("tincSnapshot"), config.get<ValueType>("DT")) == 0) {
                        wavefields->write(config.getAndCatch("snapType", 0), config.getAndCatch<std::string>("WavefieldFileName", "") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber) + ".sourceReflect", tStep, *derivatives, *modelPerShot, config.get<IndexType>("FileFormat"));
                    }
                }
                solver->resetCPML();
            }
            
            // check wavefield and seismogram for NaNs or infinite values
            if ((commShot->any(!wavefields->isFinite(dist)) || commShot->any(!receivers.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write model and break
                model->write("model_crash", config.get<IndexType>("FileFormat"));
            COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield, output model as model_crash.FILE_EXTENSION!");
            }
            if (useSourceEncode == 0) {
                sourceEst.applyOffsetMute(config, shotIndTrue, receivers);
            } else {
                sourceEst.applyOffsetMuteEncode(commShot, shotNumber, config, sourceSettingsEncode, receivers);
            }
            if (misfitType.compare("l3") == 0 || multiMisfitType.find('3') != std::string::npos) { 
                if (useSourceEncode == 0) {               
                    sourceEst.calcRefTraces(config, shotIndTrue, receivers, sourceSignalTaper);    
                } else {
                    receivers.decode(config, filenameSyn, shotNumber, sourceSettingsEncode, 0);
                    sourceEst.calcRefTracesEncode(commShot, shotNumber, config, modelCoordinates, ctx, dist, sourceSettingsEncode, receivers, sourceSignalTaper);
                }                 
            }
            if (config.get<IndexType>("useSeismogramTaper") > 1 && config.get<IndexType>("useSeismogramTaper") != 5) {                                                   
                seismogramTaper2D.apply(receivers.getSeismogramHandler()); 
            }
            seismogramTaper1D.apply(receivers.getSeismogramHandler());
            
            /* Normalize synthetic data */
            receivers.getSeismogramHandler().normalize(config.get<IndexType>("normalizeTraces"));
                            
            receivers.decode(config, filenameSyn, shotNumber, sourceSettingsEncode, 1); // for StepLengthSearch
            receivers.encode(config, filenameSyn, shotNumber, sourceSettingsEncode, 0);
            receivers.writeReceiverMark(config, shotNumber, workflow.workflowStage + 1, workflow.iteration);
            receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), filenameSyn + ".shot_" + std::to_string(shotNumber), modelCoordinates);
            
            if (useSourceEncode == 0) {
                HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Calculate misfit and adjoint sources\n");
                /* Calculate misfit of one shot */
                misfitPerIt.setValue(shotIndTrue, dataMisfit->calc(receivers, receiversTrue, shotIndTrue));
                /* Calculate adjoint sources */
                dataMisfit->calcAdjointSources(adjointSources, receivers, receiversTrue, shotIndTrue);
            } else {
                HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Calculate encode misfit and adjoint sources\n");
                /* Calculate misfit and write adjoint sources */
                dataMisfit->calcMisfitAndAdjointSources(commShot, misfitPerIt, adjointSources, receivers, receiversTrue, shotIndTrue, shotNumber, config, modelCoordinates, ctx, dist, sourceSettingsEncode, model->getVmin(), seedtime);
            }
            if (config.getAndCatch("writeAdjointSource", false))
                adjointSources.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), filenameObs + ".adjointSource" + ".shot_" + std::to_string(shotNumber), modelCoordinates);
            
            /* Calculate gradient */
            end_t_shot = common::Walltime::get();
            if (gradientKernelPerIt == 2 && decomposition == 0) {
                HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Start reflection backward in " << end_t_shot - start_t_shot << " sec.\n");
            } else {
                HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Start backward in " << end_t_shot - start_t_shot << " sec.\n");
            }
            
            gradientCalculation.run(commAll, *solver, *derivatives, receivers, sources, adjointSources, *modelPerShot, *gradientPerShot, wavefieldrecord, config, modelCoordinates, shotNumber, shotIndTrue, workflow, wavefieldTaper2D, wavefieldrecordReflect, *dataMisfit, energyPrecond, energyPrecondReflect, sourceSettingsEncode);
            
            if (!useStreamConfig) {
                *gradient += *gradientPerShot;
            } else {
                gradient->sumGradientPerShot(*model, *gradientPerShot, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotIndPerShot));
            }

            end_t_shot = common::Walltime::get();
            HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Finished in " << end_t_shot - start_t_shot << " sec.\n");

        } //end of loop over shots
        
        if (uniqueShotNos.size() == sourceSettings.size() && uniqueShotNos.size() > 1) {
            if (config.get<bool>("writeInvertedSource") && (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1)) {
                sources.getSeismogramHandler().sumShotDomain(commInterShot);
                sources.getSeismogramHandler().assignCOP();
                if (commInterShot->getRank() == 0) {
                    sources.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("sourceSeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1), modelCoordinates);
                }
                start_t_shot = common::Walltime::get();
                HOST_PRINT(commAll, "Finished sources sumShotDomain in " << start_t_shot - end_t_shot << " sec.\n", "");
            }        
            if (receivers.getNumTracesGlobal() == numShotPerSuperShot) {
                receivers.getSeismogramHandler().sumShotDomain(commInterShot);
                receivers.getSeismogramHandler().assignCOP();
                if (commInterShot->getRank() == 0) {
                    receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration), modelCoordinates);
                }
                if (config.getAndCatch("writeAdjointSource", false)) {
                    adjointSources.getSeismogramHandler().sumShotDomain(commInterShot);
                    adjointSources.getSeismogramHandler().assignCOP();
                    if (commInterShot->getRank() == 0) {
                        adjointSources.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration) + ".adjointSource", modelCoordinates);
                    }
                }
                if ((useSourceEncode == 0 && (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1)) || useSourceEncode != 0) {
                    receiversTrue.getSeismogramHandler().sumShotDomain(commInterShot, 1);
                    receiversTrue.getSeismogramHandler().assignCOP();
                    if (commInterShot->getRank() == 0) {
                        if (useSourceEncode == 0 && (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1)) {
                            receiversTrue.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1), modelCoordinates);
                        } else if (useSourceEncode != 0) {
                            receiversTrue.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration), modelCoordinates);
                        }
                    }
                    if (config.get<IndexType>("normalizeTraces") == 3 || misfitType.compare("l6") == 0 || multiMisfitType.find('6') != std::string::npos) {
                        if (commInterShot->getRank() == 0)
                            receiversTrue.getSeismogramHandler().write(5, config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1), modelCoordinates);
                    }            
                }
                if (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1) {
                    receiversStart.getSeismogramHandler().sumShotDomain(commInterShot);
                    receiversStart.getSeismogramHandler().assignCOP();
                    if (commInterShot->getRank() == 0) {
                        receiversStart.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1), modelCoordinates);
                    }
                }
                end_t_shot = common::Walltime::get();
                HOST_PRINT(commAll, "Finished receivers sumShotDomain in " << end_t_shot - start_t_shot << " sec.\n", "");
            }
        }
    
        commInterShot->sumArray(misfitPerIt.getLocalValues());
        dataMisfit->sumShotDomain(commInterShot);
        dataMisfit->addToStorage(misfitPerIt);
        misfitPerIt = 0;
        gradient->sumShotDomain(commInterShot); 
        gradient->smooth(commAll, config);

        HOST_PRINT(commAll, "\n======== Finished loop over shots " << equationType << " " << equationInd << " =========");
        HOST_PRINT(commAll, "\n=================================================\n");
        
        if (config.get<bool>("useGradientTaper"))
            gradientTaper1D.apply(*gradient);
                    
        if (inversionType == 2 || config.getAndCatch("saveCrossGradientMisfit", 0)) {
            // joint inversion with cross gradient constraint
            HOST_PRINT(commAll, "\n===========================================");
            HOST_PRINT(commAll, "\n========= calcCrossGradient " << equationType << " " << equationInd << " =========");
            HOST_PRINT(commAll, "\n===========================================\n");
            
            crossGradientDerivativeEM->calcModelDerivative(*dataMisfitEM, *modelEM, *derivativesInversionEM, configEM, modelTaper2DJoint, workflowEM);
            
            crossGradientDerivative->calcCrossGradient(*dataMisfitEM, *model, *derivativesInversionEM, configEM, modelTaper2DJoint, workflow);  
            if (config.get<IndexType>("writeGradient") == 2 && commInterShot->getRank() == 0){
                crossGradientDerivative->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".crossGradient", config.get<IndexType>("FileFormat"), workflow);
            }
                                    
            ValueType crossGradientMisfit = crossGradientDerivative->calcCrossGradientMisfit();
            dataMisfit->addToCrossGradientMisfitStorage(crossGradientMisfit);
            HOST_PRINT(commAll, "cross gradient misfit " << equationType << " " << equationInd << " = " << crossGradientMisfit << "\n"); 
            
            crossGradientDerivative->calcCrossGradientDerivative(*dataMisfitEM, *model, *derivativesInversionEM, configEM, modelTaper2DJoint, workflow);
            crossGradientDerivative->sumShotDomain(commInterShot); 
            if (config.get<IndexType>("writeGradient") == 2 && commInterShot->getRank() == 0){
                crossGradientDerivative->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".crossGradientDerivative", config.get<IndexType>("FileFormat"), workflow);
            }      
            if (inversionType == 2 && workflow.iteration > 1) {
                ValueType relativeCrossGradientMisfit = (dataMisfit->getCrossGradientMisfit(workflow.iteration - 1) - dataMisfit->getCrossGradientMisfit(workflow.iteration)) / dataMisfit->getCrossGradientMisfit(workflow.iteration - 1);
                ValueType relativeMisfit = (dataMisfit->getMisfitSum(workflow.iteration - 1) - dataMisfit->getMisfitSum(workflow.iteration)) / dataMisfit->getMisfitSum(workflow.iteration - 1);
                if (relativeMisfit > 0) {
                    weightingCrossGradient = abs(relativeCrossGradientMisfit) / (abs(relativeCrossGradientMisfit) + relativeMisfit);
                } else {
                    weightingCrossGradient = 0;
                } 
                HOST_PRINT(commAll, "weightingCrossGradient " << equationType << " " << equationInd << " = " << weightingCrossGradient << "\n");  
                HOST_PRINT(commAll, "\n===========================================\n");
                
                gradient->normalize();  
                crossGradientDerivative->normalize();  
                *crossGradientDerivative *= weightingCrossGradient;            
                *gradient += *crossGradientDerivative;
            }
        } else {
            ValueType crossGradientMisfit = 0;
            dataMisfit->addToCrossGradientMisfitStorage(crossGradientMisfit);
        }
        
        if (workflow.iteration != 0 && config.get<IndexType>("stablizingFunctionalType") != 0) {            
            // inversion with regularization constraint
            HOST_PRINT(commAll, "\n===========================================");
            HOST_PRINT(commAll, "\n=== calcStabilizingFunctionalGradient " << equationType << " " << equationInd << " ===");
            HOST_PRINT(commAll, "\n===========================================\n");
            
            stabilizingFunctionalGradient->calcStabilizingFunctionalGradient(*model, *modelPriori, config, *dataMisfit, workflow);                    
//                     stabilizingFunctionalGradient->smooth(commAll, config);  
            if (config.get<IndexType>("writeGradient") == 2 && commInterShot->getRank() == 0){
                stabilizingFunctionalGradient->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".stabilizingFunctionalGradient", config.get<IndexType>("FileFormat"), workflow);
            }  
            if ((workflow.iteration > 0) && (dataMisfit->getMisfitSum(workflow.iteration - 1) - dataMisfit->getMisfitSum(workflow.iteration) - 0.01 * dataMisfit->getMisfitSum(workflow.iteration - 1) < 0)) {
                weightingStabilizingFunctionalGradient *= 0.6;
            }
            HOST_PRINT(commAll, "weightingStabilizingFunctionalGradient = " << weightingStabilizingFunctionalGradient << "\n");  
            HOST_PRINT(commAll, "\n===========================================\n");
            
            gradient->normalize();  
            stabilizingFunctionalGradient->normalize();  
            *stabilizingFunctionalGradient *= weightingStabilizingFunctionalGradient;
            *gradient += *stabilizingFunctionalGradient;   
        }
        
        // scale function in gradientOptimization must be the final operation for gradient.
        gradientOptimization->apply(*gradient, workflow, *model, config);
        
        /* Output of gradient */
        /* only shot domain 0 writes output */
        if (config.get<IndexType>("writeGradient") > 0 && commInterShot->getRank() == 0) {
            if (useRTM == 0) {
                gradient->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1), config.get<IndexType>("FileFormat"), workflow);
            } else {
                gradient->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(maxiterations + 1), config.get<IndexType>("FileFormat"), workflow);
            }
        }
        
        SLsearch.appendToLogFile(commAll, workflow.workflowStage + 1, workflow.iteration, logFilename, dataMisfit->getMisfitSum(workflow.iteration), dataMisfit->getCrossGradientMisfit(workflow.iteration));
        dataMisfit->appendMisfitTypeShotsToFile(commAll, logFilename, workflow.workflowStage + 1, workflow.iteration);
        dataMisfit->appendMisfitPerShotToFile(commAll, logFilename, workflow.workflowStage + 1, workflow.iteration);
        dataMisfit->appendMultiMisfitsToFile(commAll, logFilename, workflow.workflowStage + 1, workflow.iteration);

        HOST_PRINT(commAll, "\nMisfit after stage " << workflow.workflowStage + 1 << ", iteration " << workflow.iteration << ": " << dataMisfit->getMisfitSum(workflow.iteration) << "\n");                
        /* --------------------------------------- */
        /* Check abort criteria for two inversions */
        /* --------------------------------------- */              
        HOST_PRINT(commAll, "\n========== Check abort criteria " << equationType << " " << equationInd << " ==========\n"); 
        breakLoop = abortCriterion.check(commAll, *dataMisfit, steplengthInit, workflow, breakLoopEM, breakLoopType);
        if (useSourceEncode != 0 || useRandomSource != 0) {
            breakLoop = false;
        }
    }
}

/*! \brief Update model of the InversionSingle
 \param commAll CommunicatorPtr
 \param dist dist
 \param model model
 \param config config
 \param modelCoordinates modelCoordinates
 \param modelCoordinatesBig modelCoordinatesBig
 \param workflow workflow
 \param dataMisfit dataMisfit
 \param crossGradientDerivative crossGradientDerivative
 \param SLsearch SLsearch
 \param modelTaper2DJoint modelTaper2DJoint
 \param maxiterations maxiterations
 \param useRTM useRTM
 \param breakLoop breakLoop
 \param ctx context
 \param seedtime seedtime
 \param inversionType inversionType
 \param equationInd equationInd
 \param breakLoopEM breakLoopEM
 \param modelEM modelEM
 \param configEM configEM
 \param modelCoordinatesEM modelCoordinatesEM
 \param workflowEM workflowEM
 \param dataMisfitEM dataMisfitEM
 \param crossGradientDerivativeEM crossGradientDerivativeEM
 \param derivativesInversionEM derivativesInversionEM
 */
template <typename ValueType>
void KITGPI::InversionSingle<ValueType>::updateModel(scai::dmemo::CommunicatorPtr commAll, scai::dmemo::DistributionPtr dist, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &model, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, KITGPI::Workflow::Workflow<ValueType> &workflow, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfit, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivative, KITGPI::StepLengthSearch<ValueType> &SLsearch, Taper::Taper2D<ValueType> modelTaper2DJoint, IndexType maxiterations, IndexType &useRTM, bool &breakLoop, scai::hmemo::ContextPtr ctx, IndexType &seedtime, IndexType inversionType, IndexType equationInd, bool &breakLoopEM, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &modelEM, KITGPI::Configuration::Configuration configEM, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesEM, KITGPI::Workflow::Workflow<ValueType> &workflowEM, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfitEM, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivativeEM, typename ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr &derivativesInversionEM)
{         
    if (inversionType != 0 && (breakLoop == false || breakLoopType == 2 || useRTM != 0)) {
        dmemo::CommunicatorPtr commShot = dist->getCommunicatorPtr();
        dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());
        SCAI_DMEMO_TASK(commShot)
        
        HOST_PRINT(commAll, "\n================================================");
        HOST_PRINT(commAll, "\n========== Start step length search " << equationType << " " << equationInd << " ======\n");
        if (config.getAndCatch("steplengthType", 2) == 1) {
            std::vector<bool> invertForParameters = workflow.getInvertForParameters();
            SLsearch.init();
            for (unsigned i=0; i<invertForParameters.size()-1; i++) {
                if (invertForParameters[i]) {
                    std::vector<bool> invertForParametersTemp(invertForParameters.size(), false);
                    invertForParametersTemp[i] = invertForParameters[i];
                    gradient->setInvertForParameters(invertForParametersTemp);
                    
                    SLsearch.run(commAll, *solver, *derivatives, sources, receivers, receiversTrue, *model, dist, config, modelCoordinates, *gradient, steplengthInit, *dataMisfit, workflow, freqFilter, sourceEst, sourceSignalTaper);

                    *gradient *= SLsearch.getSteplength();
                }
            }
            gradient->setInvertForParameters(invertForParameters);
        } else if (config.getAndCatch("steplengthType", 2) == 0 || config.getAndCatch("steplengthType", 2) == 2) {                        
            SLsearch.run(commAll, *solver, *derivatives, sources, receivers, receiversTrue, *model, dist, config, modelCoordinates, *gradient, steplengthInit, *dataMisfit, workflow, freqFilter, sourceEst, sourceSignalTaper);

            *gradient *= SLsearch.getSteplength();
        }
        HOST_PRINT(commAll, "================= Update Model " << equationType << " " << equationInd << " ============\n\n");
        /* Apply model update */
        *model -= *gradient;
        
        if (config.get<bool>("useModelThresholds"))
            model->applyThresholds(config);

        if (model->getParameterisation() == 2 || model->getParameterisation() == 1) {
            HOST_PRINT(commAll, "\n======= calcWaveModulusFromPetrophysics " << equationType << " " << equationInd << " =======\n");  
            model->calcWaveModulusFromPetrophysics();  
            if (config.get<bool>("useModelThresholds"))
                model->applyThresholds(config); 
        } else if (inversionType == 3) {
            HOST_PRINT(commAll, "\n======= calcPetrophysicsFromWaveModulus " << equationType << " " << equationInd << " =======\n");  
            model->calcPetrophysicsFromWaveModulus();
            if (config.get<bool>("useModelThresholds"))
                model->applyThresholds(config); 
            if (exchangeStrategy != 5 && exchangeStrategy != 6) {
                // case 0,1,2,3,4: self-constraint of the petrophysical relationship    
                HOST_PRINT(commAll, "\n======= calcWaveModulusFromPetrophysics " << equationType << " " << equationInd << " =======\n");  
                model->calcWaveModulusFromPetrophysics(); 
                if (config.get<bool>("useModelThresholds"))
                    model->applyThresholds(config); 
            }                    
        }
        if (commInterShot->getRank() == 0) {
            /* only shot domain 0 writes output */
            model->write((config.get<std::string>("ModelFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1)), config.get<IndexType>("FileFormat"));
        }
        
        steplengthInit *= 0.98; 
        
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "\nFinished iteration " << workflow.iteration + 1 << " in " << end_t - start_t << " sec.\n\n");
    }
}

/*! \brief Run one extra forward modelling to ensure complete and consistent output
 \param commAll CommunicatorPtr
 \param dist dist
 \param model model
 \param config config
 \param modelCoordinates modelCoordinates
 \param modelCoordinatesBig modelCoordinatesBig
 \param workflow workflow
 \param dataMisfit dataMisfit
 \param crossGradientDerivative crossGradientDerivative
 \param SLsearch SLsearch
 \param modelTaper2DJoint modelTaper2DJoint
 \param maxiterations maxiterations
 \param useRTM useRTM
 \param breakLoop breakLoop
 \param ctx context
 \param seedtime seedtime
 \param inversionType inversionType
 \param equationInd equationInd
 \param breakLoopEM breakLoopEM
 \param modelEM modelEM
 \param configEM configEM
 \param modelCoordinatesEM modelCoordinatesEM
 \param workflowEM workflowEM
 \param dataMisfitEM dataMisfitEM
 \param crossGradientDerivativeEM crossGradientDerivativeEM
 \param derivativesInversionEM derivativesInversionEM
 */
template <typename ValueType>
void KITGPI::InversionSingle<ValueType>::runExtra(scai::dmemo::CommunicatorPtr commAll, scai::dmemo::DistributionPtr dist, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &model, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, KITGPI::Workflow::Workflow<ValueType> &workflow, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfit, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivative, KITGPI::StepLengthSearch<ValueType> &SLsearch, Taper::Taper2D<ValueType> modelTaper2DJoint, IndexType maxiterations, IndexType &useRTM, bool &breakLoop, scai::hmemo::ContextPtr ctx, IndexType &seedtime, IndexType inversionType, IndexType equationInd, bool &breakLoopEM, typename Modelparameter::Modelparameter<ValueType>::ModelparameterPtr &modelEM, KITGPI::Configuration::Configuration configEM, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesEM, KITGPI::Workflow::Workflow<ValueType> &workflowEM, typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr &dataMisfitEM, typename Gradient::Gradient<ValueType>::GradientPtr &crossGradientDerivativeEM, typename ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr &derivativesInversionEM)
{          
    /* -------------------------------------------------------------------- */
    /* One extra forward modelling to ensure complete and consistent output */
    /* -------------------------------------------------------------------- */
    if (inversionType != 0 && (breakLoop == false || breakLoopType == 2 || useRTM != 0) && workflow.iteration == maxiterations - 1) {
        HOST_PRINT(commAll, "\n================ Maximum number of iterations reached " << equationType << " " << equationInd << " ================\n");
        HOST_PRINT(commAll, "== Do one more forward modelling to calculate misfit and save seismograms ==\n\n");
    
        IndexType shotNumber;
        IndexType shotIndTrue = 0; 
        IndexType shotIndIncr = 0;  
        
        dmemo::CommunicatorPtr commShot = dist->getCommunicatorPtr();
        dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());
        SCAI_DMEMO_TASK(commShot)
        
        if (!useStreamConfig) {     
            *modelPerShot = *model;
            modelPerShot->prepareForModelling(modelCoordinates, ctx, dist, commShot);
            solver->prepareForModelling(*modelPerShot, config.get<ValueType>("DT"));
        }
                
        std::vector<IndexType> uniqueShotInds = sources.getUniqueShotInds();
        std::vector<IndexType> shotIndsIncr = sources.getShotIndsIncr();
        Acquisition::writeRandomShotNosToFile(commAll, logFilename, uniqueShotNos, uniqueShotInds, workflow.workflowStage + 1, workflow.iteration + 1, useRandomSource);
        sources.writeSourceFC(commAll, config, workflow.workflowStage + 1, workflow.iteration + 1);
        sources.writeSourceEncode(commAll, config, workflow.workflowStage + 1, workflow.iteration + 1);
        dataMisfit->init(config, misfitTypeHistory, numshots, useRTM, model->getVmin(), seedtime); // in case of that random misfit function is used
        
        for (IndexType shotInd = shotDist->lb(); shotInd < shotDist->ub(); shotInd++) {
            shotIndTrue = uniqueShotInds[shotInd];
            shotIndIncr = shotIndsIncr[shotInd]; // it is not compatible with useSourceEncode != 0
            
            std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
            if (useSourceEncode == 0) {
                shotNumber = uniqueShotNos[shotIndTrue];
                Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
            } else {
                shotNumber = uniqueShotNosEncode[shotIndTrue];
                Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettingsEncode, shotNumber);
            }                    
            sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);
            if (uniqueShotNos.size() == sourceSettings.size() && uniqueShotNos.size() > 1) {
                sources.getSeismogramHandler().setShotInd(shotIndTrue, shotIndIncr);
            }

            if (useStreamConfig) {
                HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Switch to model subset\n");
                IndexType shotIndPerShot = shotIndTrue;
                if (useSourceEncode == 3) {
                    Acquisition::getuniqueShotInd(shotIndPerShot, sourceSettingsEncode, shotNumber);
                }
                model->getModelPerShot(*modelPerShot, dist, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotIndPerShot)); 
                modelPerShot->prepareForModelling(modelCoordinates, ctx, dist, commShot); 
                solver->prepareForModelling(*modelPerShot, config.get<ValueType>("DT"));
            }

            /* Read field data (or pseudo-observed data, respectively) */
            if (config.get<IndexType>("useReceiversPerShot") != 0) {
                receivers.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncode);
                receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncode);
            }
            if (uniqueShotNos.size() == sourceSettings.size() && uniqueShotNos.size() > 1 && receivers.getNumTracesGlobal() == numShotPerSuperShot) {
                receivers.getSeismogramHandler().setShotInd(shotIndTrue, shotIndIncr);
                receiversTrue.getSeismogramHandler().setShotInd(shotIndTrue, shotIndIncr);
            }
            
            if (useSourceEncode == 0) {
                receiversTrue.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".shot_" + std::to_string(shotNumber), 1);
            } else {
                receiversTrue.encode(config, config.get<std::string>("fieldSeisName"), shotNumber, sourceSettingsEncode, 1);
            }

            HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Additional forward run with " << tStepEnd << " time steps\n");

            if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0) {
                sources.getSeismogramHandler().filter(freqFilter);
                receiversTrue.getSeismogramHandler().filter(freqFilter);
            }
            
            std::string filenameSyn = config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1);
            std::string filenameObs = config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1);
            
            if (workflow.getMinOffset() != 0 || workflow.getMaxOffset() != 0 || misfitType.compare("l3") == 0 || multiMisfitType.find('3') != std::string::npos || misfitType.compare("l4") == 0 || multiMisfitType.find('4') != std::string::npos){
                if (useSourceEncode == 0) {
                    sourceEst.calcOffsetMutes(sources, receiversTrue, workflow.getMinOffset(), workflow.getMaxOffset(), shotIndTrue, modelCoordinates);
                    sourceEst.applyOffsetMute(config, shotIndTrue, receiversTrue);
                } else {
                    sourceEst.calcOffsetMutesEncode(commShot, shotNumber, config, workflow, modelCoordinates, ctx, dist, sourceSettingsEncode, receiversTrue);
                    receiversTrue.decode(config, filenameObs, shotNumber, sourceSettingsEncode, 0);
                    sourceEst.applyOffsetMuteEncode(commShot, shotNumber, config, sourceSettingsEncode, receiversTrue);
                }
                if (config.get<IndexType>("useSourceSignalInversion") == 2 || misfitType.compare("l3") == 0 || multiMisfitType.find('3') != std::string::npos) {
                    if (config.get<IndexType>("useSourceSignalTaper") == 2) {
                        sourceSignalTaper.calcCosineTaper(sources.getSeismogramHandler(), workflow.getLowerCornerFreq(), workflow.getUpperCornerFreq(), config.get<ValueType>("DT"), ctx);
                    }
                    if (useSourceEncode == 0) {
                        sourceEst.calcRefTraces(config, shotIndTrue, receiversTrue, sourceSignalTaper);
                    } else {
                        sourceEst.calcRefTracesEncode(commShot, shotNumber, config, modelCoordinates, ctx, dist, sourceSettingsEncode, receiversTrue, sourceSignalTaper);
                    }
                    sourceEst.setRefTracesToSource(sources, receiversTrue, sourceSettingsEncode, shotIndTrue, shotNumber);
                }
            }

            if (config.get<IndexType>("useSourceSignalInversion") != 0) {
                if (useSourceEncode == 0) {
                    sourceEst.applyFilter(sources, shotNumber, sourceSettings);
                } else {
                    sourceEst.applyFilter(sources, shotNumber, sourceSettingsEncode);
                }
                if (config.get<IndexType>("useSourceSignalTaper") != 0) {
                    if (config.get<IndexType>("useSourceSignalTaper") == 2) {
                        sourceSignalTaper.calcCosineTaper(sources.getSeismogramHandler(), workflow.getLowerCornerFreq(), workflow.getUpperCornerFreq(), config.get<ValueType>("DT"), ctx);
                    }
                    sourceSignalTaper.apply(sources.getSeismogramHandler());
                }
            }

            if (config.get<IndexType>("useSeismogramTaper") > 1 && config.get<IndexType>("useSeismogramTaper") != 5) {             
                seismogramTaper2D.init(receiversTrue.getSeismogramHandler());
                if (config.get<IndexType>("useSeismogramTaper") == 4) {
                    seismogramTaper2D.read(config.get<std::string>("seismogramTaperName") + ".misfitCalc.shot_" + std::to_string(shotNumber) + ".mtx");
                } else {
                    seismogramTaper2D.read(config.get<std::string>("seismogramTaperName") + ".shot_" + std::to_string(shotNumber) + ".mtx");
                }                                      
                seismogramTaper2D.apply(receiversTrue.getSeismogramHandler()); 
            }
            seismogramTaper1D.apply(receiversTrue.getSeismogramHandler());

            start_t_shot = common::Walltime::get();
            wavefields->resetWavefields();

            for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                solver->run(receivers, sources, *modelPerShot, *wavefields, *derivatives, tStep);
            }
            solver->resetCPML();

            // check wavefield and seismogram for NaNs or infinite values
            if ((commShot->any(!wavefields->isFinite(dist)) || commShot->any(!receivers.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write model and break
                model->write("model_crash", config.get<IndexType>("FileFormat"));
                COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield, output model as model_crash.FILE_EXTENSION!");
            }
            if (useSourceEncode == 0) {
                sourceEst.applyOffsetMute(config, shotIndTrue, receivers);
            } else {
                receivers.decode(config, filenameSyn, shotNumber, sourceSettingsEncode, 0);
                sourceEst.applyOffsetMuteEncode(commShot, shotNumber, config, sourceSettingsEncode, receivers);
            }
            if (misfitType.compare("l3") == 0 || multiMisfitType.find('3') != std::string::npos) {
                if (useSourceEncode == 0) {               
                    sourceEst.calcRefTraces(config, shotIndTrue, receivers, sourceSignalTaper);    
                } else {
                    sourceEst.calcRefTracesEncode(commShot, shotNumber, config, modelCoordinates, ctx, dist, sourceSettingsEncode, receivers, sourceSignalTaper);
                }   
            }
            
            if (config.get<IndexType>("useSeismogramTaper") > 1 && config.get<IndexType>("useSeismogramTaper") != 5) {                                                   
                seismogramTaper2D.apply(receivers.getSeismogramHandler()); 
            }
            seismogramTaper1D.apply(receivers.getSeismogramHandler());
            
            /* Normalize observed and synthetic data */
            if (config.get<IndexType>("normalizeTraces") == 3 || misfitType.compare("l6") == 0 || multiMisfitType.find('6') != std::string::npos) {
                // to read inverseAGC matrix.
                receiversTrue.getSeismogramHandler().read(5, config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber));         
                receivers.getSeismogramHandler().setInverseAGC(receiversTrue.getSeismogramHandler());
            }
            receivers.getSeismogramHandler().normalize(config.get<IndexType>("normalizeTraces"));
            receiversTrue.getSeismogramHandler().normalize(config.get<IndexType>("normalizeTraces"));
                
            if (config.getAndCatch("writeAdjointSource", false)) {
                receivers.decode(config, config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1), shotNumber, sourceSettingsEncode, 1);
            } else {
                receivers.decode(config, config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1), shotNumber, sourceSettingsEncode, 0);
            }
            receivers.encode(config, config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1), shotNumber, sourceSettingsEncode, 0);
            receivers.writeReceiverMark(config, shotNumber, workflow.workflowStage + 1, workflow.iteration + 1);
            receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
            
            if (useSourceEncode == 0 && shotHistory[shotIndTrue] == 1) {
                receiversTrue.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
            } else if (useSourceEncode != 0) {
                receiversTrue.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".shot_" + std::to_string(shotNumber), modelCoordinates);
            }
        
            if (useSourceEncode == 0) {
                HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Calculate misfit\n");
                /* Calculate misfit of one shot */
                misfitPerIt.setValue(shotIndTrue, dataMisfit->calc(receivers, receiversTrue, shotIndTrue));
            } else {
                HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Calculate encode misfit\n");
                /* Calculate misfit and write adjoint sources */
                IndexType seedtimeTemp = 0;
                dataMisfit->calcMisfitAndAdjointSources(commShot, misfitPerIt, adjointSources, receivers, receiversTrue, shotIndTrue, shotNumber, config, modelCoordinates, ctx, dist, sourceSettingsEncode, model->getVmin(), seedtimeTemp);
            }
            end_t_shot = common::Walltime::get();
            HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << "): Finished additional forward run\n");

        } //end of loop over shots     
        
        if (uniqueShotNos.size() == sourceSettings.size() && uniqueShotNos.size() > 1) {
            if (config.get<bool>("writeInvertedSource") && (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1)) {
                sources.getSeismogramHandler().sumShotDomain(commInterShot);
                sources.getSeismogramHandler().assignCOP();
                if (commInterShot->getRank() == 0) {
                    sources.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("sourceSeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1), modelCoordinates);
                }
                start_t_shot = common::Walltime::get();
                HOST_PRINT(commAll, "Finished sources sumShotDomain in " << start_t_shot - end_t_shot << " sec.\n", "");
            }        
            if (receivers.getNumTracesGlobal() == numShotPerSuperShot) {
                receivers.getSeismogramHandler().sumShotDomain(commInterShot);
                receivers.getSeismogramHandler().assignCOP();
                if (commInterShot->getRank() == 0) {
                    receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1), modelCoordinates);
                }
                if ((useSourceEncode == 0 && (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1)) || useSourceEncode != 0) {
                    receiversTrue.getSeismogramHandler().sumShotDomain(commInterShot, 1);
                    receiversTrue.getSeismogramHandler().assignCOP();
                    if (commInterShot->getRank() == 0) {
                        if (useSourceEncode == 0 && (workflow.iteration == 0 || shotHistory[shotIndTrue] == 1)) {
                            receiversTrue.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1), modelCoordinates);
                        } else if (useSourceEncode != 0) {
                            receiversTrue.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1), modelCoordinates);
                        }
                    }
                }
                end_t_shot = common::Walltime::get();
                HOST_PRINT(commAll, "Finished receivers sumShotDomain in " << end_t_shot - start_t_shot << " sec.\n", "");
            }
        }
            
        commInterShot->sumArray(misfitPerIt.getLocalValues());          
        dataMisfit->sumShotDomain(commInterShot);  
        dataMisfit->addToStorage(misfitPerIt);

        if (inversionType == 2 || config.getAndCatch("saveCrossGradientMisfit", 0)) {
            // joint inversion with cross gradient constraint
            HOST_PRINT(commAll, "\n===========================================");
            HOST_PRINT(commAll, "\n========= calcCrossGradient " << equationType << " " << equationInd << " =========");
            HOST_PRINT(commAll, "\n===========================================\n");
            
            crossGradientDerivativeEM->calcModelDerivative(*dataMisfitEM, *modelEM, *derivativesInversionEM, configEM, modelTaper2DJoint, workflowEM);
            
            crossGradientDerivative->calcCrossGradient(*dataMisfitEM, *model, *derivativesInversionEM, configEM, modelTaper2DJoint, workflow);  
            if (config.get<IndexType>("writeGradient") == 2 && commInterShot->getRank() == 0){
                crossGradientDerivative->write(gradname + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration + 1) + ".crossGradient", config.get<IndexType>("FileFormat"), workflow);
            }
                                    
            ValueType crossGradientMisfit = crossGradientDerivative->calcCrossGradientMisfit();
            dataMisfit->addToCrossGradientMisfitStorage(crossGradientMisfit);
            HOST_PRINT(commAll, "cross gradient misfit " << equationType << " " << equationInd << " = " << crossGradientMisfit << "\n");
        } else {
            ValueType crossGradientMisfit = 0;
            dataMisfit->addToCrossGradientMisfitStorage(crossGradientMisfit);
        }
        
        SLsearch.appendToLogFile(commAll, workflow.workflowStage + 1, workflow.iteration + 1, logFilename, dataMisfit->getMisfitSum(workflow.iteration + 1), dataMisfit->getCrossGradientMisfit(workflow.iteration + 1));
        dataMisfit->appendMisfitTypeShotsToFile(commAll, logFilename, workflow.workflowStage + 1, workflow.iteration + 1);
        dataMisfit->appendMisfitPerShotToFile(commAll, logFilename, workflow.workflowStage + 1, workflow.iteration + 1);
        dataMisfit->appendMultiMisfitsToFile(commAll, logFilename, workflow.workflowStage + 1, workflow.iteration + 1);

        HOST_PRINT(commAll, "\nMisfit after stage " << workflow.workflowStage + 1 << ", iteration " << workflow.iteration + 1 << ": " << dataMisfit->getMisfitSum(workflow.iteration + 1) << "\n");  
        
        if (breakLoopType == 0 && breakLoopEM == true && (exchangeStrategy == 3 || exchangeStrategy == 5)) {
            breakLoop = true;
            breakLoopEM = false;
            workflow.iteration = 0;
        }                 
        if (gradientKernel == 4) {
            useRTM = 1;
            workflow.iteration--;
        }     
    } // end extra Seismic forward modelling
}

template class KITGPI::InversionSingle<double>;
template class KITGPI::InversionSingle<float>;
