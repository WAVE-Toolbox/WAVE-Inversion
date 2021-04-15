
/*! \brief Calculate misfit 
 *
 *
 \param solverEM Forward solver
 \param derivativesEM Derivatives matrices
 \param receiversEM Receivers
 \param sourcesEM Sources 
 \param modelEM Model for the finite-difference simulation
 \param wavefieldsEM Wavefields
 \param configEM Configuration
 \param scaledGradientEM Misfit gradient 
 \param steplength Steplength
 */
template <typename ValueType>
ValueType KITGPI::StepLengthSearch<ValueType>::calcMisfit(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolverEM<ValueType> &solverEM, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivativesEM, KITGPI::Acquisition::ReceiversEM<ValueType> &receiversEM, std::vector<Acquisition::sourceSettings<ValueType>> &sourceSettingsEM, KITGPI::Acquisition::ReceiversEM<ValueType> &receiversTrueEM, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Wavefields::WavefieldsEM<ValueType> &wavefieldsEM, KITGPI::Configuration::Configuration configEM, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesEM, KITGPI::Gradient::GradientEM<ValueType> &scaledGradientEM, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, ValueType steplength, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM, KITGPI::Filter::Filter<ValueType> const &freqFilterEM, KITGPI::SourceEstimation<ValueType> const &sourceEstEM, KITGPI::Taper::Taper1D<ValueType> const &sourceSignalTaperEM, std::vector<scai::IndexType> uniqueShotNosRandEM)
{
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    std::string equationTypeEM = configEM.get<std::string>("equationType");
    std::transform(equationTypeEM.begin(), equationTypeEM.end(), equationTypeEM.begin(), ::tolower); 
    scai::dmemo::DistributionPtr distEM;
    if(equationTypeEM.compare("tmem") == 0 || equationTypeEM.compare("viscotmem") == 0){
        distEM = wavefieldsEM.getRefEZ().getDistributionPtr();
    } else {
        distEM = wavefieldsEM.getRefEX().getDistributionPtr();      
    }
    scai::dmemo::CommunicatorPtr commShot = distEM->getCommunicatorPtr();
    scai::dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    SCAI_DMEMO_TASK(commShot)

    Acquisition::SourcesEM<ValueType> sourcesEM;
    std::vector<scai::IndexType> uniqueShotNosEM;
    Acquisition::calcuniqueShotNo(uniqueShotNosEM, sourceSettingsEM);
    IndexType numshotsEM = uniqueShotNosEM.size();
    std::shared_ptr<const dmemo::BlockDistribution> shotDistEM;
    if (configEM.get<IndexType>("useRandSource") != 0) {  
        IndexType numShotDomainsEM = configEM.get<IndexType>("NumShotDomains");
        shotDistEM = dmemo::blockDistribution(numShotDomainsEM, commInterShot);
    } else {
        shotDistEM = dmemo::blockDistribution(numshotsEM, commInterShot);
    }
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    //     double start_t, end_t; /* For timing */
    IndexType tStepEndEM = static_cast<IndexType>((configEM.get<ValueType>("T") / configEM.get<ValueType>("DT")) + 0.5);
    int testShotIncr = configEM.get<int>("testShotIncr");

    scai::lama::DenseVector<ValueType> misfitTestEM(numshotsEM, 0, ctx);

    // Implement a (virtual) copy constructor in the abstract base class to simplify the following code -> virtual constructor idiom!
    typename KITGPI::Modelparameter::ModelparameterEM<ValueType>::ModelparameterPtr testmodelEM(KITGPI::Modelparameter::FactoryEM<ValueType>::Create(equationTypeEM));
    *testmodelEM = modelEM;
    testmodelEM->prepareForInversion(configEM, commShot);
    
    typename KITGPI::Modelparameter::ModelparameterEM<ValueType>::ModelparameterPtr testmodelPerShotEM(KITGPI::Modelparameter::FactoryEM<ValueType>::Create(equationTypeEM));
    testmodelPerShotEM->init(configEM, ctx, distEM, modelCoordinatesEM);  // init is necessary for modelPerShot to allocate distribution.
    bool useStreamConfigEM = configEM.get<bool>("useStreamConfig");
    Acquisition::Coordinates<ValueType> modelCoordinatesBigEM;
    std::vector<Acquisition::coordinate3D> cutCoordinatesEM;
    if (useStreamConfigEM) {
        KITGPI::Configuration::Configuration configBigEM(configEM.get<std::string>("streamConfigFilename"));
        modelCoordinatesBigEM.init(configBigEM);
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBigEM;
        sourcesEM.getAcquisitionSettings(configBigEM, sourceSettingsBigEM);
        Acquisition::getCutCoord(cutCoordinatesEM, sourceSettingsBigEM);
    }
    
    typename KITGPI::Gradient::GradientEM<ValueType>::GradientPtr testgradient(KITGPI::Gradient::FactoryEM<ValueType>::Create(equationTypeEM));
    *testgradient = scaledGradientEM;

    /* Update modelEM */
    *testgradient *= steplength;
    *testmodelEM -= *testgradient;
    
    if (configEM.get<bool>("useModelThresholds"))
        testmodelEM->applyThresholds(configEM);

    if (testmodelEM->getParameterisation() == 2 || testmodelEM->getParameterisation() == 1) {
        HOST_PRINT(commAll, "\n======= calcWaveModulusFromPetrophysics test =====\n");  
        testmodelEM->calcWaveModulusFromPetrophysics();  
    }
    
    if (configEM.get<bool>("useModelThresholds"))
        testmodelEM->applyThresholds(configEM); 

    if (!useStreamConfigEM) {
        testmodelEM->prepareForModelling(modelCoordinatesEM, ctx, distEM, commShot);
        solverEM.prepareForModelling(*testmodelEM, configEM.get<ValueType>("DT"));
    }

    IndexType shotNumber;  
    IndexType shotIndTrue = 0;  
    // later it should be possible to select only a subset of shots for the step length search
    for (IndexType shotInd = shotDistEM->lb(); shotInd < shotDistEM->ub(); shotInd += testShotIncr) {
        if (configEM.get<IndexType>("useRandSource") == 0) {  
            shotNumber = uniqueShotNosEM[shotInd];
            shotIndTrue = shotInd;
        } else {
            shotNumber = uniqueShotNosRandEM[shotInd];
            Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNosEM, shotNumber);
        }

        if (useStreamConfigEM) {
            HOST_PRINT(commShot, "Switch to model shot: " << shotIndTrue + 1 << " of " << numshotsEM << "\n");
            testmodelEM->getModelPerShot(*testmodelPerShotEM, modelCoordinatesEM, modelCoordinatesBigEM, cutCoordinatesEM.at(shotIndTrue));    
            testmodelPerShotEM->prepareForModelling(modelCoordinatesEM, ctx, distEM, commShot); 
            solverEM.prepareForModelling(*testmodelPerShotEM, configEM.get<ValueType>("DT"));
        }

        wavefieldsEM.resetWavefields();

        HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshotsEM << ": Start Test Forward\n");

        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
        Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettingsEM, shotNumber);
        sourcesEM.init(sourceSettingsShot, configEM, modelCoordinatesEM, ctx, distEM);

        if (!useStreamConfigEM) {
            CheckParameter::checkNumericalArtifactsAndInstabilities<ValueType>(configEM, sourceSettingsShot, *testmodelEM, modelCoordinatesEM, shotNumber);
        } else {
            CheckParameter::checkNumericalArtifactsAndInstabilities<ValueType>(configEM, sourceSettingsShot, *testmodelPerShotEM, modelCoordinatesEM, shotNumber);
        }
        
        if (configEM.get<bool>("useReceiversPerShot")) {
            receiversEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber);
            receiversTrueEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber);
        } else if (configEM.get<bool>("useReceiverMark")) {
            receiversEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber, numshotsEM);
            receiversTrueEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber, numshotsEM);
        }

        receiversTrueEM.getSeismogramHandler().read(configEM.get<IndexType>("SeismogramFormat"), configEM.get<std::string>("FieldSeisName") + ".shot_" + std::to_string(shotNumber), 1);

        if (workflowEM.getLowerCornerFreq() != 0.0 || workflowEM.getUpperCornerFreq() != 0.0) {
            sourcesEM.getSeismogramHandler().filter(freqFilterEM);
            receiversTrueEM.getSeismogramHandler().filter(freqFilterEM);
        }
        
        Taper::Taper2D<ValueType> seismogramTaper2DEM;
        Taper::Taper1D<ValueType> seismogramTaper1DEM;
        seismogramTaper1DEM.init(std::make_shared<dmemo::NoDistribution>(tStepEndEM), ctx, 1);
        seismogramTaper1DEM.calcTimeDampingTaper(workflowEM.getTimeDampingFactor(), configEM.get<ValueType>("DT"));  
        if (configEM.get<IndexType>("useSeismogramTaper") == 2) {
            seismogramTaper2DEM.init(receiversTrueEM.getSeismogramHandler());
            seismogramTaper2DEM.read(configEM.get<std::string>("seismogramTaperName") + ".shot_" + std::to_string(shotNumber) + ".mtx");                       
            seismogramTaper2DEM.apply(receiversTrueEM.getSeismogramHandler());  
        }
        seismogramTaper1DEM.apply(receiversTrueEM.getSeismogramHandler());

        if (configEM.get<bool>("useSourceSignalInversion")) {
            sourceEstEM.applyFilter(sourcesEM, shotIndTrue);
            if (configEM.get<bool>("useSourceSignalTaper"))
                sourceSignalTaperEM.apply(sourcesEM.getSeismogramHandler());
        }
        
        if (!useStreamConfigEM) {
            for (IndexType tStep = 0; tStep < tStepEndEM; tStep++) {
                solverEM.run(receiversEM, sourcesEM, *testmodelEM, wavefieldsEM, derivativesEM, tStep);
            }
        } else {
            for (IndexType tStep = 0; tStep < tStepEndEM; tStep++) {
                solverEM.run(receiversEM, sourcesEM, *testmodelPerShotEM, wavefieldsEM, derivativesEM, tStep);
            }
        }

        // check wavefield and seismogram for NaNs or infinite values
        if ((commShot->any(!wavefieldsEM.isFinite(distEM)) || commShot->any(!receiversEM.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write modelEM and break
            testmodelEM->write("model_crash", configEM.get<IndexType>("FileFormat"));
            COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield for modelEM in steplength search, output modelEM as model_crash.FILE_EXTENSION!");
        }

        if (configEM.get<IndexType>("useSeismogramTaper") == 2) {                                                   
            seismogramTaper2DEM.apply(receiversEM.getSeismogramHandler()); 
        }
        seismogramTaper1DEM.apply(receiversEM.getSeismogramHandler());

        /* Normalize observed and synthetic data */
        if (configEM.get<bool>("normalizeTraces")) {
            receiversEM.getSeismogramHandler().normalize();
            receiversTrueEM.getSeismogramHandler().normalize();
        }

        misfitTestEM.setValue(shotIndTrue, dataMisfitEM.calc(receiversEM, receiversTrueEM, shotIndTrue));
    }

    commInterShot->sumArray(misfitTestEM.getLocalValues());

    HOST_PRINT(commAll, "\n======== Finished loop over test shots =========");
    HOST_PRINT(commAll, "\n================================================\n");

    return misfitTestEM.sum();
}
