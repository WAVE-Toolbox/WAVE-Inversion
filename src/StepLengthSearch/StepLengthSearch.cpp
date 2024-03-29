#include "StepLengthSearch.hpp"
#include <iomanip>
#include <string>
#include <IO/IO.hpp>

using scai::IndexType;

/*! \brief Find the optimal steplength
 *
 *
 \param solver Forward solver
 \param derivatives Derivatives matrices
 \param receivers Receivers
 \param sources Sources 
 \param model Model for the finite-difference simulation
 \param dist Distribution
 \param config Configuration
 \param scaledGradient Misfit gradient 
 \param steplengthInit Initial steplength
 \param currentMisfit Current misfit
 */
template <typename ValueType>
void KITGPI::StepLengthSearch<ValueType>::run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Sources<ValueType> sources, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, ValueType steplengthInit, KITGPI::Misfit::Misfit<ValueType> &currentMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Filter::Filter<ValueType> const &freqFilter, KITGPI::SourceEstimation<ValueType> sourceEst, KITGPI::Taper::Taper1D<ValueType> sourceSignalTaper)
{
    scaledGradient.printInvertForParameters(commAll);
    if (steplengthType == 0) {
        HOST_PRINT(commAll, "Constant steplength: " << steplengthInit << " \n");
        steplengthOptimum = steplengthInit;
    } else if (steplengthType == 1) {
        this->runLineSearch(commAll, solver, derivatives, sources, receivers, receiversTrue, model, dist, config, modelCoordinates, scaledGradient, steplengthInit, currentMisfit, workflow, freqFilter, sourceEst, sourceSignalTaper);
    } else if (steplengthType == 2) {
        this->runParabolicSearch(commAll, solver, derivatives, sources, receivers, receiversTrue, model, dist, config, modelCoordinates, scaledGradient, steplengthInit, currentMisfit, workflow, freqFilter, sourceEst, sourceSignalTaper);
    }
}

/*! \brief Find the optimal steplength
 *
 *
 \param solver Forward solver
 \param derivatives Derivatives matrices
 \param receivers Receivers
 \param sources Sources 
 \param model Model for the finite-difference simulation
 \param dist Distribution
 \param config Configuration
 \param scaledGradient Misfit gradient 
 \param steplengthInit Initial steplength
 \param currentMisfit Current misfit
 */
template <typename ValueType>
void KITGPI::StepLengthSearch<ValueType>::runLineSearch(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Sources<ValueType> sources, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, ValueType steplengthInit, KITGPI::Misfit::Misfit<ValueType> &currentMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Filter::Filter<ValueType> const &freqFilter, KITGPI::SourceEstimation<ValueType> sourceEst, KITGPI::Taper::Taper1D<ValueType> sourceSignalTaper)
{
    double start_t, end_t; /* For timing */

    scai::dmemo::CommunicatorPtr commShot = dist->getCommunicatorPtr();
    scai::dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings = sources.getSourceSettings();
    std::vector<scai::IndexType> uniqueShotNos;
    Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
    IndexType numshots = 1;
    IndexType numShotDomains = config.get<IndexType>("NumShotDomains");
    Common::checkNumShotDomains(numShotDomains, commAll);
    IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    if (useSourceEncode == 0) {
        numshots = uniqueShotNos.size();
    } else {
        numshots = numShotDomains;
    }
    std::shared_ptr<const dmemo::BlockDistribution> shotDist;
    if (config.getAndCatch("useRandomSource", 0) != 0) {  
        shotDist = dmemo::blockDistribution(numShotDomains, commInterShot);
    } else {
        shotDist = dmemo::blockDistribution(numshots, commInterShot);
    }
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);   
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower); 

    scai::IndexType numRelaxationMechanisms = config.get<IndexType>("numRelaxationMechanisms");
    typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType));
    wavefields->init(ctx, dist, numRelaxationMechanisms);

    typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(KITGPI::Misfit::Factory<ValueType>::Create(config.get<std::string>("misfitType")));
    *dataMisfit = currentMisfit;

    ValueType misfitTestSum;

    /* ------------------------------------------- */
    /* Set values for step length search           */
    /* ------------------------------------------- */
    ValueType scalingFactor = config.get<ValueType>("scalingFactor");
    ValueType steplengthMin = config.get<ValueType>("steplengthMin");
    ValueType steplengthMax = config.get<ValueType>("steplengthMax");

    /* ------------------------------------------------------------------------------------------------------ */
    /* See details in Meles et al., 2010 or Pica et al., 1990                                                 */
    /* ------------------------------------------------------------------------------------------------------ */

    start_t = scai::common::Walltime::get();

    HOST_PRINT(commAll, "\nEstimation of steplength by line search\n");
    misfitTestSum = this->calcMisfit(commAll, solver, derivatives, sources, receivers, receiversTrue, model, *wavefields, config, modelCoordinates, scaledGradient, *dataMisfit, steplengthInit, workflow, freqFilter, sourceEst, sourceSignalTaper);
    HOST_PRINT(commAll, "\nOptimum step length in line search: " << steplengthOptimum << "\n");

    steplengthGuess = steplengthInit;
    /* Check if accepted step length is smaller than maximally allowed step length */
    if (steplengthOptimum > steplengthMax) {
        steplengthOptimum = steplengthMax;
        HOST_PRINT(commAll, "\nVariable steplengthMax used to update the model");
    } else if (steplengthOptimum < steplengthMin) {
        steplengthOptimum = steplengthMin;
        HOST_PRINT(commAll, "\nVariable steplengthMin used to update the model");
    }
    if (steplengthOptimum > steplengthGuess * scalingFactor) {
        steplengthOptimum = steplengthGuess * scalingFactor;
        HOST_PRINT(commAll, "\nVariable steplengthInit and scalingFactor used to update the model");
    }
    if (std::isnan(steplengthOptimum))
        steplengthOptimum = steplengthMin;
    std::vector<bool> invertForParameters = scaledGradient.getInvertForParameters();
    for (unsigned i=0; i<invertForParameters.size(); i++) {
        if (invertForParameters[i]) {
            steplengthLine.setValue(i, steplengthOptimum);
            misfitLine.setValue(i, misfitTestSum);
            
            /* Output of the test step lengths and the corresponding misfit */
            HOST_PRINT(commAll, "\nSteplength " << i + 1 << ": " << steplengthGuess << ", Corresponding misfit: " << misfitLine.getValue(i));
            break;
        }
    }
    HOST_PRINT(commAll, "\nOptimum step length: " << steplengthOptimum << "\n");

    end_t = scai::common::Walltime::get();
    HOST_PRINT(commAll, "\nFinished step length search in " << end_t - start_t << " sec.\n\n\n");
}

/*! \brief Find the optimal steplength
 *
 *
 \param solver Forward solver
 \param derivatives Derivatives matrices
 \param receivers Receivers
 \param sources Sources 
 \param model Model for the finite-difference simulation
 \param dist Distribution
 \param config Configuration
 \param scaledGradient Misfit gradient 
 \param steplengthInit Initial steplength
 \param currentMisfit Current misfit
 */
template <typename ValueType>
void KITGPI::StepLengthSearch<ValueType>::runParabolicSearch(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Sources<ValueType> sources, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, ValueType steplengthInit, KITGPI::Misfit::Misfit<ValueType> &currentMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Filter::Filter<ValueType> const &freqFilter, KITGPI::SourceEstimation<ValueType> sourceEst, KITGPI::Taper::Taper1D<ValueType> sourceSignalTaper)
{
    double start_t, end_t; /* For timing */

    scai::dmemo::CommunicatorPtr commShot = dist->getCommunicatorPtr();
    scai::dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings = sources.getSourceSettings();
    std::vector<scai::IndexType> uniqueShotNos;
    Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
    IndexType numshots = 1;
    IndexType numShotDomains = config.get<IndexType>("NumShotDomains");
    Common::checkNumShotDomains(numShotDomains, commAll);
    IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    if (useSourceEncode == 0) {
        numshots = uniqueShotNos.size();
    } else {
        numshots = numShotDomains;
    }
    std::shared_ptr<const dmemo::BlockDistribution> shotDist;
    if (config.getAndCatch("useRandomSource", 0) != 0) {  
        shotDist = dmemo::blockDistribution(numShotDomains, commInterShot);
    } else {
        shotDist = dmemo::blockDistribution(numshots, commInterShot);
    }
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);   
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower); 
    int testShotIncr = config.get<int>("testShotIncr");

    scai::IndexType numRelaxationMechanisms = config.get<IndexType>("numRelaxationMechanisms");
    typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType));
    wavefields->init(ctx, dist, numRelaxationMechanisms);

    typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(KITGPI::Misfit::Factory<ValueType>::Create(config.get<std::string>("misfitType")));
    *dataMisfit = currentMisfit;

    ValueType misfitTestSum;
    ValueType steplength;

    /* ------------------------------------------- */
    /* Set values for step length search           */
    /* ------------------------------------------- */
    stepCalcCount = 1;                                // number of forward calculations to find a proper (steplength, misfit) pair
    int maxStepCalc = config.get<int>("maxStepCalc"); // maximum number of forward calculations to find a proper (steplength, misfit) pair
    ValueType scalingFactor = config.get<ValueType>("scalingFactor");
    ValueType steplengthMin = config.get<ValueType>("steplengthMin");
    ValueType steplengthMax = config.get<ValueType>("steplengthMax");

    /* Save three pairs (steplength, misfit) for the parabolic fit */
    steplengthParabola.setSameValue(3, 0);
    misfitParabola.setSameValue(3, 0);

    /* ------------------------------------------------------------------------------------------------------ */
    /* Try to find a second step length which gives a smaller misfit and a third step length which gives again a larger misfit.
       The first step length is zero so the misfit from the gradient calculation can be used to save computation time.
       
       Separate 4 cases: 
       1) misfit 2 < misfit 1 AND misfit 3 > misfit 2: SL = minimum of parabola or SL = SLmax
       2) misfit 2 < misfit 1 AND misfit 3 < misfit 2: SL = SL3 or SL = SLmax 
       3) misfit 2 > misfit 1 AND misfit 3 < misfit 2 AND misfit 3 < misfit 1: SL = SL3 
       4) misfit 2 > misfit 1 AND misfit 3 > misfit 1: SL = very small SL */
    /* ------------------------------------------------------------------------------------------------------ */

    start_t = scai::common::Walltime::get();

    step2ok = false; // true if second step length decreases the misfit
    step3ok = false; // true if second step length decreases the misfit AND third step length increases the misfit relative to second step length */

    /* --- Save first step length (steplength=0 is used to save computational time) --- */
    std::vector<IndexType> uniqueShotInds = sources.getUniqueShotInds();
    misfitTestSum = 0; 
    IndexType shotIndTrue = 0;  
    for (IndexType shotInd = shotDist->lb(); shotInd < shotDist->ub(); shotInd += testShotIncr) {
        shotIndTrue = uniqueShotInds[shotInd];
        misfitTestSum += currentMisfit.getMisfitIt(workflow.iteration).getValue(shotIndTrue);
    }
    misfitTestSum = commInterShot->sum(misfitTestSum);

    misfitParabola.setValue(0, misfitTestSum);

    /* --- Save second step length (initial step length) in any case --- */
    HOST_PRINT(commAll, "\nEstimation of 2nd steplength, forward test run no. " << stepCalcCount << " of maximum " << maxStepCalc << "\n");
    misfitTestSum = this->calcMisfit(commAll, solver, derivatives, sources, receivers, receiversTrue, model, *wavefields, config, modelCoordinates, scaledGradient, *dataMisfit, steplengthInit, workflow, freqFilter, sourceEst, sourceSignalTaper);
    steplengthParabola.setValue(1, steplengthInit);
    misfitParabola.setValue(1, misfitTestSum);
    if (misfitParabola.getValue(0) > misfitParabola.getValue(1)) {
        step2ok = true;
    }

    /* --- Search for a third step length - case: misfit was DECREASED Case 1/2 --- */
    steplength = steplengthInit;
    while (step2ok == true && stepCalcCount < maxStepCalc) {
        HOST_PRINT(commAll, "\nEstimation of 3rd steplength, forward test run no. " << stepCalcCount + 1 << " of maximum " << maxStepCalc << "\n");
        steplength *= scalingFactor;
        misfitTestSum = this->calcMisfit(commAll, solver, derivatives, sources, receivers, receiversTrue, model, *wavefields, config, modelCoordinates, scaledGradient, *dataMisfit, steplength, workflow, freqFilter, sourceEst, sourceSignalTaper);
        steplengthParabola.setValue(2, steplength);
        misfitParabola.setValue(2, misfitTestSum);
        if (misfitTestSum > misfitParabola.getValue(1)) {
            step3ok = true;
            stepCalcCount += 1;
            break;
        } else {
            stepCalcCount += 1;
        }
    }

    /* --- Search for a third step length - case: misfit was INCREASED Case 3/4 --- */
    steplength = steplengthInit;
    while (step2ok == false && stepCalcCount < maxStepCalc) {
        HOST_PRINT(commAll, "Estimation of 3rd steplength, forward test run no. " << stepCalcCount + 1 << " of maximum " << maxStepCalc << "\n");
        steplength /= scalingFactor;
        misfitTestSum = this->calcMisfit(commAll, solver, derivatives, sources, receivers, receiversTrue, model, *wavefields, config, modelCoordinates, scaledGradient, *dataMisfit, steplength, workflow, freqFilter, sourceEst, sourceSignalTaper);
        steplengthParabola.setValue(2, steplength);
        misfitParabola.setValue(2, misfitTestSum);
        if (misfitTestSum < misfitParabola.getValue(0)) {
            step3ok = true;
            stepCalcCount += 1;
            break;
        } else {
            stepCalcCount += 1;
        }
    }

    /* Output of the three test step lengths and the corresponding misfit */
    for (int i = 0; i < 3; i++) {
        HOST_PRINT(commAll, "\nSteplength " << i + 1 << ": " << steplengthParabola.getValue(i) << ", Corresponding misfit: " << misfitParabola.getValue(i));
    }

    /* Set optimum step length */
    if (step2ok == true && step3ok == true) {
        steplengthOptimum = parabolicFit(steplengthParabola, misfitParabola);
        HOST_PRINT(commAll, "\nApply parabolic fit");
    } else if (step2ok == true && step3ok == false) {
        if (maxStepCalc==1){
            steplengthOptimum = steplengthParabola.getValue(1);}
        else {
            steplengthOptimum = steplengthParabola.getValue(2);}
        HOST_PRINT(commAll, "\nVariable steplengthInit and scalingFactor used to update the model");
    } else if (step2ok == false && step3ok == true) {
        steplengthOptimum = parabolicFit(steplengthParabola, misfitParabola); //case 3: parabolic fit instead of steplengthParabola.getValue(2)
        HOST_PRINT(commAll, "\nApply parabolic fit");
    } else if (step2ok == false && step3ok == false) {
        steplengthOptimum = steplengthMin;
        HOST_PRINT(commAll, "\nVariable steplengthMin used to update the model");
    }

    /* Check if accepted step length is smaller than maximally allowed step length */
    if (steplengthOptimum > steplengthMax) {
        steplengthOptimum = steplengthMax;
        HOST_PRINT(commAll, "\nVariable steplengthMax used to update the model");
    } else if (steplengthOptimum < steplengthMin) {
        steplengthOptimum = steplengthMin;
        HOST_PRINT(commAll, "\nVariable steplengthMin used to update the model");
    }

    if (std::isnan(steplengthOptimum))
        steplengthOptimum = steplengthMin;

    HOST_PRINT(commAll, "\nOptimum step length: " << steplengthOptimum << "\n");

    end_t = scai::common::Walltime::get();
    HOST_PRINT(commAll, "\nFinished step length search in " << end_t - start_t << " sec.\n\n\n");
}

/*! \brief Parabolic fit
 *
 *
 \param xValues Vector of three values on the x-axis
 \param xValues Vector of three values on the y-axis
 */
template <typename ValueType>
ValueType KITGPI::StepLengthSearch<ValueType>::parabolicFit(scai::lama::DenseVector<ValueType> const &xValues, scai::lama::DenseVector<ValueType> const &yValues)
{
    SCAI_ASSERT(xValues.size() == 3, "xvalues must contain 3 values!");
    SCAI_ASSERT(yValues.size() == 3, "yvalues must contain 3 values!");
    ValueType steplengthExtremum;
    scai::lama::DenseMatrix<ValueType> A(3, 3);
    scai::lama::DenseMatrix<ValueType> invA;
    scai::lama::DenseVector<ValueType> coeff; // does the size need to be specified?
    scai::lama::DenseVector<ValueType> vectorOnes(3, 1, scai::hmemo::Context::getContextPtr());
    scai::lama::DenseVector<ValueType> xValuesPow2 = xValues;
    xValuesPow2 *= xValues;

    A.setColumn(xValuesPow2, 0, scai::common::BinaryOp::COPY);
    A.setColumn(xValues, 1, scai::common::BinaryOp::COPY);
    A.setColumn(vectorOnes, 2, scai::common::BinaryOp::COPY);
    invA.invert(A);

    coeff = invA * yValues;
    steplengthExtremum = -coeff.getValue(1) / (2 * coeff.getValue(0));

    return steplengthExtremum;
}

/*! \brief Calculate misfit 
 *
 *
 \param solver Forward solver
 \param derivatives Derivatives matrices
 \param receivers Receivers
 \param sources Sources 
 \param model Model for the finite-difference simulation
 \param wavefields Wavefields
 \param config Configuration
 \param scaledGradient Misfit gradient 
 \param steplength Steplength
 */
template <typename ValueType>
ValueType KITGPI::StepLengthSearch<ValueType>::calcMisfit(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Sources<ValueType> sources, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, ValueType steplength, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Filter::Filter<ValueType> const &freqFilter, KITGPI::SourceEstimation<ValueType> sourceEst, KITGPI::Taper::Taper1D<ValueType> sourceSignalTaper)
{
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
    std::string equationType = config.get<std::string>("equationType");
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower); 
    bool isSeismic = Common::checkEquationType<ValueType>(equationType);   
    scai::dmemo::DistributionPtr dist;
    if (isSeismic) {
        if(equationType.compare("sh") == 0 || equationType.compare("viscosh") == 0){
            dist = wavefields.getRefVZ().getDistributionPtr();
        } else {
            dist = wavefields.getRefVX().getDistributionPtr();      
        }
    } else {
        if(equationType.compare("tmem") == 0 || equationType.compare("viscotmem") == 0){
            dist = wavefields.getRefEZ().getDistributionPtr();
        } else {
            dist = wavefields.getRefEX().getDistributionPtr();      
        }
    }
    scai::dmemo::CommunicatorPtr commShot = dist->getCommunicatorPtr();
    scai::dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    SCAI_DMEMO_TASK(commShot)

    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings = sources.getSourceSettings();
    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsEncode;
    std::vector<scai::IndexType> uniqueShotNos;
    std::vector<scai::IndexType> uniqueShotNosEncode;
    IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    IndexType numshots = 1;
    IndexType numShotDomains = config.get<IndexType>("NumShotDomains");
    Common::checkNumShotDomains(numShotDomains, commAll);
    Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
    if (useSourceEncode == 0) {
        numshots = uniqueShotNos.size();
    } else {
        sourceSettingsEncode = sources.getSourceSettingsEncode();
        Acquisition::calcuniqueShotNo(uniqueShotNosEncode, sourceSettingsEncode);
        numshots = numShotDomains;
    }
    bool useStreamConfig = config.getAndCatch("useStreamConfig", false);
    Acquisition::Coordinates<ValueType> modelCoordinatesBig;
    std::vector<Acquisition::coordinate3D> cutCoordinates;
    if (useStreamConfig) {
        KITGPI::Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
        modelCoordinatesBig.init(configBig);
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;
        ValueType shotIncr = config.getAndCatch("shotIncr", 0.0);
        sources.getAcquisitionSettings(config, shotIncr);
        sourceSettingsBig = sources.getSourceSettings(); 
        Acquisition::getCutCoord(config, cutCoordinates, sourceSettingsBig, modelCoordinates, modelCoordinatesBig);
    }
    
    std::shared_ptr<const dmemo::BlockDistribution> shotDist;
    if (config.getAndCatch("useRandomSource", 0) != 0) {  
        shotDist = dmemo::blockDistribution(numShotDomains, commInterShot);
    } else {
        shotDist = dmemo::blockDistribution(numshots, commInterShot);
    }
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    Acquisition::Receivers<ValueType> receiversLast;
    if (config.get<IndexType>("useReceiversPerShot") == 0) {
        receiversLast.init(config, modelCoordinates, ctx, dist);
    }
    IndexType numshotsIncr = sourceSettings.size();
    IndexType numShotPerSuperShot = ceil(ValueType(numshotsIncr) / numShotDomains); 
    if (uniqueShotNos.size() == sourceSettings.size() && uniqueShotNos.size() > 1 && receivers.getNumTracesGlobal() == numShotPerSuperShot) {
        receiversLast.getSeismogramHandler().allocateCOP(numshots, tStepEnd);
    }

    int testShotIncr = config.get<int>("testShotIncr");

    scai::lama::DenseVector<ValueType> misfitTest(numshots, 0, ctx);

    // Implement a (virtual) copy constructor in the abstract base class to simplify the following code -> virtual constructor idiom!
    typename KITGPI::Modelparameter::Modelparameter<ValueType>::ModelparameterPtr testmodel(KITGPI::Modelparameter::Factory<ValueType>::Create(equationType));
    *testmodel = model;
    testmodel->prepareForInversion(config, commShot);
    
    typename KITGPI::Modelparameter::Modelparameter<ValueType>::ModelparameterPtr testmodelPerShot(KITGPI::Modelparameter::Factory<ValueType>::Create(equationType));
    
    typename KITGPI::Gradient::Gradient<ValueType>::GradientPtr testgradient(KITGPI::Gradient::Factory<ValueType>::Create(equationType));
    *testgradient = scaledGradient;
    std::vector<bool> invertForParameters = scaledGradient.getInvertForParameters();
    testgradient->setInvertForParameters(invertForParameters);

    /* Update model */
    *testgradient *= steplength;
    *testmodel -= *testgradient;
    
    if (config.get<bool>("useModelThresholds"))
        testmodel->applyThresholds(config);

    if (testmodel->getParameterisation() == 2 || testmodel->getParameterisation() == 1) {
        HOST_PRINT(commAll, "\n======= calcWaveModulusFromPetrophysics test =====\n");  
        testmodel->calcWaveModulusFromPetrophysics();  
    }
    
    if (config.get<bool>("useModelThresholds"))
        testmodel->applyThresholds(config); 

    if (!useStreamConfig) {
        testmodel->prepareForModelling(modelCoordinates, ctx, dist, commShot);
        solver.prepareForModelling(*testmodel, config.get<ValueType>("DT"));
    }

    std::vector<IndexType> uniqueShotInds = sources.getUniqueShotInds();
    std::vector<IndexType> shotIndsIncr = sources.getShotIndsIncr();
    IndexType shotNumber;  
    IndexType shotIndTrue = 0;  
    IndexType shotIndIncr = 0;  
    ValueType numerator = 0;
    ValueType denominator = 0;
    // later it should be possible to select only a subset of shots for the step length search
    for (IndexType shotInd = shotDist->lb(); shotInd < shotDist->ub(); shotInd += testShotIncr) {
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

        if (!useStreamConfig) {
            CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(config, sourceSettingsShot, *testmodel, modelCoordinates, shotNumber);
        } else {
            HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << commInterShot->getRank() << ", index " << shotIndTrue + 1 << " of " << numshots << "): Switch to model subset\n");
            IndexType shotIndPerShot = shotIndTrue;
            if (useSourceEncode == 3) {
                Acquisition::getuniqueShotInd(shotIndPerShot, sourceSettingsEncode, shotNumber);
            }
            testmodel->getModelPerShot(*testmodelPerShot, dist, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotIndPerShot));    
            testmodelPerShot->prepareForModelling(modelCoordinates, ctx, dist, commShot); 
            solver.prepareForModelling(*testmodelPerShot, config.get<ValueType>("DT"));
            
            CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(config, sourceSettingsShot, *testmodelPerShot, modelCoordinates, shotNumber);
        }

        HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << commInterShot->getRank() << ", index " << shotIndTrue + 1 << " of " << numshots << "): Start Test Forward\n");
        
        if (config.get<IndexType>("useReceiversPerShot") != 0) {
            receivers.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncode);
            receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncode);
            receiversLast.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncode);
        }
        if (uniqueShotNos.size() == sourceSettings.size() && uniqueShotNos.size() > 1 && receivers.getNumTracesGlobal() == numShotPerSuperShot) {
            receiversLast.getSeismogramHandler().setShotInd(shotIndTrue, shotIndIncr);
        }

        if (useSourceEncode == 0) {
            receiversTrue.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("fieldSeisName") + ".shot_" + std::to_string(shotNumber), 1);
        } else {
            receiversTrue.encode(config, config.get<std::string>("fieldSeisName"), shotNumber, sourceSettingsEncode, 1);
        }
        if (useSourceEncode == 0 || receivers.getNumTracesGlobal() == numShotPerSuperShot) {
            receiversLast.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration) + ".shot_" + std::to_string(shotNumber));
        } else {
            receiversLast.encode(config, config.get<std::string>("SeismogramFilename") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".It_" + std::to_string(workflow.iteration), shotNumber, sourceSettingsEncode, 1);
        }

        if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0) {
            sources.getSeismogramHandler().filter(freqFilter);
            receiversTrue.getSeismogramHandler().filter(freqFilter);
        }
                
        if (useSourceEncode == 0) {
            sourceEst.applyOffsetMute(config, shotIndTrue, receiversTrue);
            sourceEst.applyOffsetMute(config, shotIndTrue, receiversLast);
        } else {
            receiversTrue.decode(config, "", shotNumber, sourceSettingsEncode, 0);
            sourceEst.applyOffsetMuteEncode(commShot, shotNumber, config, sourceSettingsEncode, receiversTrue);            
            sourceEst.applyOffsetMuteEncode(commShot, shotNumber, config, sourceSettingsEncode, receiversLast);
        }
        if (config.get<IndexType>("useSourceSignalInversion") == 2 || dataMisfit.getMisfitTypeShots().getValue(shotIndTrue) == 3) {
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
        
        Taper::Taper2D<ValueType> seismogramTaper2D;
        Taper::Taper1D<ValueType> seismogramTaper1D;
        seismogramTaper1D.init(std::make_shared<dmemo::NoDistribution>(tStepEnd), ctx, 1);
        seismogramTaper1D.calcTimeDampingTaper(workflow.getTimeDampingFactor(), config.get<ValueType>("DT"));  
        if (config.get<IndexType>("useSeismogramTaper") > 1) {
            seismogramTaper2D.init(receiversTrue.getSeismogramHandler());
            if (config.get<IndexType>("useSeismogramTaper") == 5) {
                seismogramTaper2D.calcCosineTaper(receiversTrue.getSeismogramHandler(), workflow.getUpperCornerFreq(), workflow.getUpperCornerFreq(), config, shotIndTrue, ctx);
            } else if (config.get<IndexType>("useSeismogramTaper") == 4) {
                seismogramTaper2D.read(config.get<std::string>("seismogramTaperName") + ".misfitCalc.shot_" + std::to_string(shotNumber) + ".mtx");
            } else {
                seismogramTaper2D.read(config.get<std::string>("seismogramTaperName") + ".shot_" + std::to_string(shotNumber) + ".mtx");
            }
            seismogramTaper2D.apply(receiversTrue.getSeismogramHandler());  
        }
        seismogramTaper1D.apply(receiversTrue.getSeismogramHandler());

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

        wavefields.resetWavefields();
        
        if (!useStreamConfig) {
            for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                solver.run(receivers, sources, *testmodel, wavefields, derivatives, tStep);
            }
        } else {
            for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                solver.run(receivers, sources, *testmodelPerShot, wavefields, derivatives, tStep);
            }
        }
        solver.resetCPML();

        // check wavefield and seismogram for NaNs or infinite values
        if ((commShot->any(!wavefields.isFinite(dist)) || commShot->any(!receivers.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write model and break
            testmodel->write("model_crash", config.get<IndexType>("FileFormat"));
            COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield for model in steplength search, output model as model_crash.FILE_EXTENSION!");
        }
        if (useSourceEncode == 0) {
            sourceEst.applyOffsetMute(config, shotIndTrue, receivers);
        } else {
            receivers.decode(config, "", shotNumber, sourceSettingsEncode, 0);
            sourceEst.applyOffsetMuteEncode(commShot, shotNumber, config, sourceSettingsEncode, receivers);
        }
        if (dataMisfit.getMisfitTypeShots().getValue(shotIndTrue) == 3) {
            if (useSourceEncode == 0) {
                sourceEst.calcRefTraces(config, shotIndTrue, receivers, sourceSignalTaper);
            } else {
                sourceEst.calcRefTracesEncode(commShot, shotNumber, config, modelCoordinates, ctx, dist, sourceSettingsEncode, receivers, sourceSignalTaper);
            }
        }

        if (config.get<IndexType>("useSeismogramTaper") > 1) {                                                   
            seismogramTaper2D.apply(receivers.getSeismogramHandler()); 
        }
        seismogramTaper1D.apply(receivers.getSeismogramHandler());

        /* Normalize observed and synthetic data */
        if (config.get<IndexType>("normalizeTraces") == 3 || dataMisfit.getMisfitTypeShots().getValue(shotIndTrue) == 6) {
            // to read inverseAGC matrix.
            receiversTrue.getSeismogramHandler().read(5, config.get<std::string>("fieldSeisName") + ".stage_" + std::to_string(workflow.workflowStage + 1) + ".shot_" + std::to_string(shotNumber));      
            receivers.getSeismogramHandler().setInverseAGC(receiversTrue.getSeismogramHandler());      
        }
        receivers.getSeismogramHandler().normalize(config.get<IndexType>("normalizeTraces"));
        receiversTrue.getSeismogramHandler().normalize(config.get<IndexType>("normalizeTraces"));
        receivers.decode(config, "", shotNumber, sourceSettingsEncode, 0);
        receivers.encode(config, "", shotNumber, sourceSettingsEncode, 0);

        if (useSourceEncode == 0) {
            /* Calculate misfit of one shot */
            misfitTest.setValue(shotIndTrue, dataMisfit.calc(receivers, receiversTrue, shotIndTrue));
        } else {
            /* Calculate misfit and write adjoint sources */
            KITGPI::Acquisition::Receivers<ValueType> adjointSources;
            IndexType seedtimeTemp;
            dataMisfit.calcMisfitAndAdjointSources(commShot, misfitTest, adjointSources, receivers, receiversTrue, shotIndTrue, shotNumber, config, modelCoordinates, ctx, dist, sourceSettingsEncode, testmodel->getVmin(), seedtimeTemp);
        }
        
        if (steplengthType == 1) {    
            KITGPI::Acquisition::Seismogram<ValueType> seismogramSyn;
            KITGPI::Acquisition::Seismogram<ValueType> seismogramSynLast;
            KITGPI::Acquisition::Seismogram<ValueType> seismogramObs;    
            KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerSyn;
            KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerSynLast;
            KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerObs;
            seismoHandlerSyn = receivers.getSeismogramHandler();
            seismoHandlerSynLast = receiversLast.getSeismogramHandler();
            seismoHandlerObs = receiversTrue.getSeismogramHandler();
            for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
                seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
                seismogramSynLast = seismoHandlerSynLast.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
                seismogramObs = seismoHandlerObs.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
                seismogramObs -= seismogramSynLast;
                seismogramSyn -= seismogramSynLast;
                denominator += seismogramSyn.getData().l2Norm() * seismogramSyn.getData().l2Norm();
                seismogramObs *= seismogramSyn;
                scai::lama::DenseVector<ValueType> traceSum = seismogramObs.getTraceSum();
                numerator += traceSum.sum();
            }  
        }
    }
    if (steplengthType == 1) {    
        steplengthOptimum = steplength * numerator / denominator;
    }

    commInterShot->sumArray(misfitTest.getLocalValues());

    HOST_PRINT(commAll, "\n======== Finished loop over test shots =========");
    HOST_PRINT(commAll, "\n================================================\n");

    return misfitTest.sum();
}

/*! \brief Initialize log-file
 *
 *
 \param comm Communicator
 \param logFilename Name of log-file
 */
template <typename ValueType>
void KITGPI::StepLengthSearch<ValueType>::initLogFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, std::string misfitType, scai::IndexType setSteplengthType, scai::IndexType setInvertNumber, scai::IndexType setSaveCrossGradientMisfit)
{
    steplengthType = setSteplengthType;
    invertNumber = setInvertNumber;
    saveCrossGradientMisfit = setSaveCrossGradientMisfit;
    int myRank = comm->getRank();
    if (myRank == MASTERGPI) {
        logFile.open(logFilename);
        logFile << "# Step length log file  \n";
        logFile << "# MisfitType = " << misfitType << "\n";
        logFile << "# Iteration 0 shows misfit of initial model of each workflow stage (only first two columns and last column is meaningful here)\n";
        if (steplengthType == 1) {
            logFile << "# Stage | Iteration";
            for (int i = 0; i < invertNumber; i++) {
                logFile << " | optimum step length " << i+1;                
            }
            logFile << " | step length guess";
            for (int i = 0; i < invertNumber; i++) {
                logFile << " | misfit of slg" << i+1;                    
            }
            if (saveCrossGradientMisfit) {
                logFile << " | final misfit of all shots | misfit of cross gradient\n";
            } else {
                logFile << " | final misfit of all shots\n";
            }
        } else if (steplengthType == 2) {
            if (saveCrossGradientMisfit) {
                logFile << "# Stage | Iteration | optimum step length | forward calculations | step length guess 1 | step length guess 2 | step length guess 3 | misfit of slg1 | misfit of slg2 | misfit of slg3 | final misfit of all shots | misfit of cross gradient\n";
            } else {
                logFile << "# Stage | Iteration | optimum step length | forward calculations | step length guess 1 | step length guess 2 | step length guess 3 | misfit of slg1 | misfit of slg2 | misfit of slg3 | final misfit of all shots\n";
            }
        }
        logFile.close();
    }
}

/*! \brief Append results of steplength search to log-file
 *
 *
 \param comm Communicator
 \param iteration Iteration count
 */
template <typename ValueType>
void KITGPI::StepLengthSearch<ValueType>::appendToLogFile(scai::dmemo::CommunicatorPtr comm, IndexType workflowStage, scai::IndexType iteration, std::string logFilename, ValueType misfitSum, ValueType crossGradientMisfit)
{
    int myRank = comm->getRank();
    if (myRank == MASTERGPI && steplengthType == 1) {
        if (iteration == 0) {
            std::string filename(logFilename);
            logFile.open(filename, std::ios_base::app);
            logFile << std::scientific;
            logFile << std::setw(5) << workflowStage << std::setw(11) << iteration;
            for (int i = 0; i < invertNumber; i++) {
                logFile << std::setw(24) << 0;                
            }
            logFile << std::setw(21) << steplengthGuess;
            for (int i = 0; i < invertNumber; i++) {
                logFile << std::setw(17) << 0;                    
            }
            if (saveCrossGradientMisfit) {
                logFile << std::setw(22) << misfitSum << std::setw(22) << crossGradientMisfit << "\n";
            } else {
                logFile << std::setw(22) << misfitSum << "\n";
            }
            logFile.close();
        } else {
            std::string filename(logFilename);
            logFile.open(filename, std::ios_base::app);
            logFile << std::scientific;
            logFile << std::setw(5) << workflowStage << std::setw(11) << iteration;
            for (int i = 0; i < invertNumber; i++) {
                ValueType steplengthLine1 = steplengthLine.getValue(i);
                logFile << std::setw(24) << steplengthLine1;                
            }
            logFile << std::setw(21) << steplengthGuess;
            for (int i = 0; i < invertNumber; i++) {
                ValueType misfitLine1 = misfitLine.getValue(i);
                logFile << std::setw(17) << misfitLine1;                    
            }
            if (saveCrossGradientMisfit) {
                logFile << std::setw(22) << misfitSum << std::setw(22) << crossGradientMisfit << "\n";
            } else {
                logFile << std::setw(22) << misfitSum << "\n";
            }
            logFile.close();
        }
    } else if (myRank == MASTERGPI && steplengthType == 2) {
        /* The following temporaries are only necessary because of a problem with LAMA: e.g. steplengthParabola.getValue(0).getValue<ValueType>() produces an error */
        ValueType steplengthParabola0 = steplengthParabola.getValue(0);
        ValueType steplengthParabola1 = steplengthParabola.getValue(1);
        ValueType steplengthParabola2 = steplengthParabola.getValue(2);
        ValueType misfitParabola0 = misfitParabola.getValue(0);
        ValueType misfitParabola1 = misfitParabola.getValue(1);
        ValueType misfitParabola2 = misfitParabola.getValue(2);

        if (iteration == 0) {
            std::string filename(logFilename);
            logFile.open(filename, std::ios_base::app);
            logFile << std::scientific;
            logFile << std::setw(5) << workflowStage << std::setw(11) << iteration << std::setw(22) << 0 << std::setw(17) << 0 << std::setw(30) << 0 << std::setw(22) << 0 << std::setw(22) << 0 << std::setw(19) << 0 << std::setw(17) << 0 << std::setw(17) << 0;
            if (saveCrossGradientMisfit) {
                logFile << std::setw(22) << misfitSum << std::setw(22) << crossGradientMisfit << "\n";
            } else {
                logFile << std::setw(22) << misfitSum << "\n";
            }
            logFile.close();
        } else {
            std::string filename(logFilename);
            logFile.open(filename, std::ios_base::app);
            logFile << std::scientific;
            logFile << std::setw(5) << workflowStage << std::setw(11) << iteration << std::setw(22) << steplengthOptimum << std::setw(17) << stepCalcCount << std::setw(30) << steplengthParabola0 << std::setw(22) << steplengthParabola1 << std::setw(22) << steplengthParabola2 << std::setw(19) << misfitParabola0 << std::setw(17) << misfitParabola1 << std::setw(17) << misfitParabola2;
            if (saveCrossGradientMisfit) {
                logFile << std::setw(22) << misfitSum << std::setw(22) << crossGradientMisfit << "\n";
            } else {
                logFile << std::setw(22) << misfitSum << "\n";
            }
            logFile.close();
        }
    }
}

/*! \brief Append the run time to log-file
 \param comm Communicator
 \param run_t Run time
 */
template <typename ValueType>
void KITGPI::StepLengthSearch<ValueType>::appendRunTimeToLogFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, double run_t)
{
    int myRank = comm->getRank();
    if (myRank == MASTERGPI) {
        std::string filename(logFilename);
        logFile.open(filename, std::ios_base::app);
        logFile << std::scientific;
        logFile << "\n# Total runtime of WAVE-Inversion in second : " << run_t;
        logFile.close();
    }
}

/*! \brief Get optimum steplength
 *
 *
 */
template <typename ValueType>
ValueType const &KITGPI::StepLengthSearch<ValueType>::getSteplength()
{
    return (steplengthOptimum);
}

/*! \brief Initialize steplength and misfit vectors
 *
 *
 */
template <typename ValueType>
void KITGPI::StepLengthSearch<ValueType>::init()
{
    if (steplengthType == 1) {
        scai::lama::DenseVector<ValueType> steplengthLineTemp(invertNumber, 0);
        scai::lama::DenseVector<ValueType> misfitLineTemp(invertNumber, 0);
        steplengthLine = steplengthLineTemp;
        misfitLine = misfitLineTemp;
    }
}

template class KITGPI::StepLengthSearch<double>;
template class KITGPI::StepLengthSearch<float>;
