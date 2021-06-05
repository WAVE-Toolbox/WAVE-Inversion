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
void KITGPI::StepLengthSearch<ValueType>::run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, std::vector<Acquisition::sourceSettings<ValueType>> &sourceSettings, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, ValueType steplengthInit, scai::lama::DenseVector<ValueType> currentMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Filter::Filter<ValueType> const &freqFilter, KITGPI::SourceEstimation<ValueType> const &sourceEst, KITGPI::Taper::Taper1D<ValueType> const &sourceSignalTaper, std::vector<scai::IndexType> uniqueShotNosRand, std::vector<std::string> misfitTypeShots)
{
    double start_t, end_t; /* For timing */

    scai::dmemo::CommunicatorPtr commShot = dist->getCommunicatorPtr();
    scai::dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    std::vector<scai::IndexType> uniqueShotNos;
    Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
    IndexType numshots = uniqueShotNos.size();
    std::shared_ptr<const dmemo::BlockDistribution> shotDist;
    if (config.getAndCatch("useRandomSource", 0) != 0) {  
        IndexType numShotDomains = config.get<IndexType>("NumShotDomains");
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

    typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType));
    wavefields->init(ctx, dist);

    typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(KITGPI::Misfit::Factory<ValueType>::Create(config.get<std::string>("misfitType")));
    dataMisfit->setMisfitTypeShots(misfitTypeShots);

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
    misfitTestSum = 0;
    IndexType shotNumber;  
    IndexType shotIndTrue = 0;  
    for (IndexType shotInd = shotDist->lb(); shotInd < shotDist->ub(); shotInd += testShotIncr) {
        if (config.getAndCatch("useRandomSource", 0) == 0) { 
            shotIndTrue = shotInd;
        } else {
            shotNumber = uniqueShotNosRand[shotInd];
            Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNos, shotNumber);
        }
        misfitTestSum += currentMisfit.getValue(shotIndTrue);
    }
    misfitTestSum = commInterShot->sum(misfitTestSum);

    misfitParabola.setValue(0, misfitTestSum);

    /* --- Save second step length (initial step length) in any case --- */
    HOST_PRINT(commAll, "\nEstimation of 2nd steplength, forward test run no. " << stepCalcCount << " of maximum " << maxStepCalc << "\n");
    misfitTestSum = this->calcMisfit(commAll, solver, derivatives, receivers, sourceSettings, receiversTrue, model, *wavefields, config, modelCoordinates, scaledGradient, *dataMisfit, steplengthInit, workflow, freqFilter, sourceEst, sourceSignalTaper, uniqueShotNosRand);
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
        misfitTestSum = this->calcMisfit(commAll, solver, derivatives, receivers, sourceSettings, receiversTrue, model, *wavefields, config, modelCoordinates, scaledGradient, *dataMisfit, steplength, workflow, freqFilter, sourceEst, sourceSignalTaper, uniqueShotNosRand);
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
        misfitTestSum = this->calcMisfit(commAll, solver, derivatives, receivers, sourceSettings, receiversTrue, model, *wavefields, config, modelCoordinates, scaledGradient, *dataMisfit, steplength, workflow, freqFilter, sourceEst, sourceSignalTaper, uniqueShotNosRand);
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
    }

    if (std::isnan(steplengthOptimum))
        steplengthOptimum = steplengthMin;

    HOST_PRINT(commAll, "\nOptimum step length: " << steplengthOptimum << "\n");

    end_t = scai::common::Walltime::get();
    HOST_PRINT(commAll, "\nFinished step length search in  " << end_t - start_t << " sec.\n\n\n");
}

/*! \brief Find the optimal steplength
 *
 *
 \param solverEM Forward solver
 \param derivativesEM Derivatives matrices
 \param receiversEM Receivers
 \param sourcesEM Sources 
 \param modelEM Model for the finite-difference simulation
 \param distEM Distribution
 \param configEM Configuration
 \param scaledGradientEM Misfit gradient 
 \param steplengthInitEM Initial steplength
 \param currentMisfitEM Current misfit
 */
template <typename ValueType>
void KITGPI::StepLengthSearch<ValueType>::run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolverEM<ValueType> &solverEM, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivativesEM, KITGPI::Acquisition::ReceiversEM<ValueType> &receiversEM, std::vector<Acquisition::sourceSettings<ValueType>> &sourceSettingsEM, KITGPI::Acquisition::ReceiversEM<ValueType> &receiversTrueEM, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, scai::dmemo::DistributionPtr distEM, KITGPI::Configuration::Configuration configEM, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesEM, KITGPI::Gradient::GradientEM<ValueType> &scaledGradientEM, ValueType steplengthInitEM, scai::lama::DenseVector<ValueType> currentMisfitEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM, KITGPI::Filter::Filter<ValueType> const &freqFilterEM, KITGPI::SourceEstimation<ValueType> const &sourceEstEM, KITGPI::Taper::Taper1D<ValueType> const &sourceSignalTaperEM, std::vector<scai::IndexType> uniqueShotNosRandEM, std::vector<std::string> misfitTypeShotsEM)
{
    double start_t, end_t; /* For timing */

    scai::dmemo::CommunicatorPtr commShot = distEM->getCommunicatorPtr();
    scai::dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    std::vector<scai::IndexType> uniqueShotNosEM;
    Acquisition::calcuniqueShotNo(uniqueShotNosEM, sourceSettingsEM);
    IndexType numshotsEM = uniqueShotNosEM.size();
    std::shared_ptr<const dmemo::BlockDistribution> shotDistEM;
    if (configEM.getAndCatch("useRandomSource", 0) != 0) {  
        IndexType numShotDomainsEM = configEM.get<IndexType>("NumShotDomains");
        shotDistEM = dmemo::blockDistribution(numShotDomainsEM, commInterShot);
    } else {
        shotDistEM = dmemo::blockDistribution(numshotsEM, commInterShot);
    }
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    std::string dimensionEM = configEM.get<std::string>("dimension");
    std::string equationTypeEM = configEM.get<std::string>("equationType");
    std::transform(dimensionEM.begin(), dimensionEM.end(), dimensionEM.begin(), ::tolower);   
    std::transform(equationTypeEM.begin(), equationTypeEM.end(), equationTypeEM.begin(), ::tolower); 
    int testShotIncr = configEM.get<int>("testShotIncr");

    typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr wavefieldsEM(KITGPI::Wavefields::FactoryEM<ValueType>::Create(dimensionEM, equationTypeEM));
    wavefieldsEM->init(ctx, distEM);

    typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr dataMisfitEM(KITGPI::Misfit::Factory<ValueType>::Create(configEM.get<std::string>("misfitType")));
    dataMisfitEM->setMisfitTypeShots(misfitTypeShotsEM);

    ValueType misfitTestSum;
    ValueType steplength;

    /* ------------------------------------------- */
    /* Set values for step length search           */
    /* ------------------------------------------- */
    stepCalcCount = 1;                                // number of forward calculations to find a proper (steplength, misfit) pair
    int maxStepCalc = configEM.get<int>("maxStepCalc"); // maximum number of forward calculations to find a proper (steplength, misfit) pair
    ValueType scalingFactor = configEM.get<ValueType>("scalingFactor");
    ValueType steplengthMin = configEM.get<ValueType>("steplengthMin");
    ValueType steplengthMax = configEM.get<ValueType>("steplengthMax");

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
    misfitTestSum = 0;
    IndexType shotNumber;  
    IndexType shotIndTrue = 0;  
    for (IndexType shotInd = shotDistEM->lb(); shotInd < shotDistEM->ub(); shotInd += testShotIncr) {
        if (configEM.getAndCatch("useRandomSource", 0) == 0) { 
            shotIndTrue = shotInd;
        } else {
            shotNumber = uniqueShotNosRandEM[shotInd];
            Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNosEM, shotNumber);
        }
        misfitTestSum += currentMisfitEM.getValue(shotIndTrue);
    }
    misfitTestSum = commInterShot->sum(misfitTestSum);

    misfitParabola.setValue(0, misfitTestSum);

    /* --- Save second step length (initial step length) in any case --- */
    HOST_PRINT(commAll, "\nEstimation of 2nd steplength, forward test run no. " << stepCalcCount << " of maximum " << maxStepCalc << "\n");
    misfitTestSum = this->calcMisfit(commAll, solverEM, derivativesEM, receiversEM, sourceSettingsEM, receiversTrueEM, modelEM, *wavefieldsEM, configEM, modelCoordinatesEM, scaledGradientEM, *dataMisfitEM, steplengthInitEM, workflowEM, freqFilterEM, sourceEstEM, sourceSignalTaperEM, uniqueShotNosRandEM);
    steplengthParabola.setValue(1, steplengthInitEM);
    misfitParabola.setValue(1, misfitTestSum);
    if (misfitParabola.getValue(0) > misfitParabola.getValue(1)) {
        step2ok = true;
    }

    /* --- Search for a third step length - case: misfit was DECREASED Case 1/2 --- */
    steplength = steplengthInitEM;
    while (step2ok == true && stepCalcCount < maxStepCalc) {
        HOST_PRINT(commAll, "\nEstimation of 3rd steplength, forward test run no. " << stepCalcCount + 1 << " of maximum " << maxStepCalc << "\n");
        steplength *= scalingFactor;
        misfitTestSum = this->calcMisfit(commAll, solverEM, derivativesEM, receiversEM, sourceSettingsEM, receiversTrueEM, modelEM, *wavefieldsEM, configEM, modelCoordinatesEM, scaledGradientEM, *dataMisfitEM, steplength, workflowEM, freqFilterEM, sourceEstEM, sourceSignalTaperEM, uniqueShotNosRandEM);
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
    steplength = steplengthInitEM;
    while (step2ok == false && stepCalcCount < maxStepCalc) {
        HOST_PRINT(commAll, "Estimation of 3rd steplength, forward test run no. " << stepCalcCount + 1 << " of maximum " << maxStepCalc << "\n");
        steplength /= scalingFactor;
        misfitTestSum = this->calcMisfit(commAll, solverEM, derivativesEM, receiversEM, sourceSettingsEM, receiversTrueEM, modelEM, *wavefieldsEM, configEM, modelCoordinatesEM, scaledGradientEM, *dataMisfitEM, steplength, workflowEM, freqFilterEM, sourceEstEM, sourceSignalTaperEM, uniqueShotNosRandEM);
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
    }

    if (std::isnan(steplengthOptimum))
        steplengthOptimum = steplengthMin;

    HOST_PRINT(commAll, "\nOptimum step length: " << steplengthOptimum << "\n");

    end_t = scai::common::Walltime::get();
    HOST_PRINT(commAll, "\nFinished step length search in  " << end_t - start_t << " sec.\n\n\n");
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
ValueType KITGPI::StepLengthSearch<ValueType>::calcMisfit(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, std::vector<Acquisition::sourceSettings<ValueType>> &sourceSettings, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, ValueType steplength, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Filter::Filter<ValueType> const &freqFilter, KITGPI::SourceEstimation<ValueType> const &sourceEst, KITGPI::Taper::Taper1D<ValueType> const &sourceSignalTaper, std::vector<scai::IndexType> uniqueShotNosRand)
{
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    std::string equationType = config.get<std::string>("equationType");
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower); 
    scai::dmemo::DistributionPtr dist;
    if(equationType.compare("sh") == 0){
        dist = wavefields.getRefVZ().getDistributionPtr();
    } else {
        dist = wavefields.getRefVX().getDistributionPtr();      
    }
    scai::dmemo::CommunicatorPtr commShot = dist->getCommunicatorPtr();
    scai::dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    SCAI_DMEMO_TASK(commShot)

    Acquisition::Sources<ValueType> sources;
    std::vector<scai::IndexType> uniqueShotNos;
    Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
    IndexType numshots = uniqueShotNos.size();
    std::shared_ptr<const dmemo::BlockDistribution> shotDist;
    if (config.getAndCatch("useRandomSource", 0) != 0) {  
        IndexType numShotDomains = config.get<IndexType>("NumShotDomains");
        shotDist = dmemo::blockDistribution(numShotDomains, commInterShot);
    } else {
        shotDist = dmemo::blockDistribution(numshots, commInterShot);
    }
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    //     double start_t, end_t; /* For timing */
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
    int testShotIncr = config.get<int>("testShotIncr");

    scai::lama::DenseVector<ValueType> misfitTest(numshots, 0, ctx);

    // Implement a (virtual) copy constructor in the abstract base class to simplify the following code -> virtual constructor idiom!
    typename KITGPI::Modelparameter::Modelparameter<ValueType>::ModelparameterPtr testmodel(KITGPI::Modelparameter::Factory<ValueType>::Create(equationType));
    *testmodel = model;
    testmodel->prepareForInversion(config, commShot);
    
    typename KITGPI::Modelparameter::Modelparameter<ValueType>::ModelparameterPtr testmodelPerShot(KITGPI::Modelparameter::Factory<ValueType>::Create(equationType));
    testmodelPerShot->init(config, ctx, dist, modelCoordinates);  // init is necessary for modelPerShot to allocate distribution.
    bool useStreamConfig = config.getAndCatch("useStreamConfig", false);
    Acquisition::Coordinates<ValueType> modelCoordinatesBig;
    std::vector<Acquisition::coordinate3D> cutCoordinates;
    if (useStreamConfig) {
        KITGPI::Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
        modelCoordinatesBig.init(configBig);
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;
        sources.getAcquisitionSettings(configBig, sourceSettingsBig);
        Acquisition::getCutCoord(cutCoordinates, sourceSettingsBig);
    }
    
    typename KITGPI::Gradient::Gradient<ValueType>::GradientPtr testgradient(KITGPI::Gradient::Factory<ValueType>::Create(equationType));
    *testgradient = scaledGradient;

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

    IndexType shotNumber;  
    IndexType shotIndTrue = 0;  
    // later it should be possible to select only a subset of shots for the step length search
    for (IndexType shotInd = shotDist->lb(); shotInd < shotDist->ub(); shotInd += testShotIncr) {
        if (config.getAndCatch("useRandomSource", 0) == 0) {  
            shotNumber = uniqueShotNos[shotInd];
            shotIndTrue = shotInd;
        } else {
            shotNumber = uniqueShotNosRand[shotInd];
            Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNos, shotNumber);
        }

        if (useStreamConfig) {
            HOST_PRINT(commShot, "Switch to model shot: " << shotIndTrue + 1 << " of " << numshots << "\n");
            testmodel->getModelPerShot(*testmodelPerShot, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotIndTrue));    
            testmodelPerShot->prepareForModelling(modelCoordinates, ctx, dist, commShot); 
            solver.prepareForModelling(*testmodelPerShot, config.get<ValueType>("DT"));
        }

        wavefields.resetWavefields();

        HOST_PRINT(commShot, "Shot " << shotIndTrue + 1 << " of " << numshots << ": Start Test Forward\n");

        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
        Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
        sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);

        if (!useStreamConfig) {
            CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(config, sourceSettingsShot, *testmodel, modelCoordinates, shotNumber);
        } else {
            CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(config, sourceSettingsShot, *testmodelPerShot, modelCoordinates, shotNumber);
        }
        
        if (config.get<IndexType>("useReceiversPerShot") == 1) {
            receivers.init(config, modelCoordinates, ctx, dist, shotNumber);
            receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber);
        } else if (config.get<IndexType>("useReceiversPerShot") == 2) {
            receivers.init(config, modelCoordinates, ctx, dist, shotNumber, numshots);
            receiversTrue.init(config, modelCoordinates, ctx, dist, shotNumber, numshots);
        }

        receiversTrue.getSeismogramHandler().read(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("FieldSeisName") + ".shot_" + std::to_string(shotNumber), 1);

        if (workflow.getLowerCornerFreq() != 0.0 || workflow.getUpperCornerFreq() != 0.0) {
            sources.getSeismogramHandler().filter(freqFilter);
            receiversTrue.getSeismogramHandler().filter(freqFilter);
        }
        
        Taper::Taper2D<ValueType> seismogramTaper2D;
        Taper::Taper1D<ValueType> seismogramTaper1D;
        seismogramTaper1D.init(std::make_shared<dmemo::NoDistribution>(tStepEnd), ctx, 1);
        seismogramTaper1D.calcTimeDampingTaper(workflow.getTimeDampingFactor(), config.get<ValueType>("DT"));  
        if (config.get<IndexType>("useSeismogramTaper") == 2 || config.get<IndexType>("useSeismogramTaper") == 3) {
            seismogramTaper2D.init(receiversTrue.getSeismogramHandler());
            seismogramTaper2D.read(config.get<std::string>("seismogramTaperName") + ".shot_" + std::to_string(shotNumber) + ".mtx");                       
            seismogramTaper2D.apply(receiversTrue.getSeismogramHandler());  
        }
        seismogramTaper1D.apply(receiversTrue.getSeismogramHandler());

        if (config.get<bool>("useSourceSignalInversion")) {
            sourceEst.applyFilter(sources, shotIndTrue);
            if (config.get<bool>("useSourceSignalTaper"))
                sourceSignalTaper.apply(sources.getSeismogramHandler());
        }
        
        if (!useStreamConfig) {
            for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                solver.run(receivers, sources, *testmodel, wavefields, derivatives, tStep);
            }
        } else {
            for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
                solver.run(receivers, sources, *testmodelPerShot, wavefields, derivatives, tStep);
            }
        }

        // check wavefield and seismogram for NaNs or infinite values
        if ((commShot->any(!wavefields.isFinite(dist)) || commShot->any(!receivers.getSeismogramHandler().isFinite())) && (commInterShot->getRank() == 0)){ // if any processor returns isfinite=false, write model and break
            testmodel->write("model_crash", config.get<IndexType>("FileFormat"));
            COMMON_THROWEXCEPTION("Infinite or NaN value in seismogram or/and velocity wavefield for model in steplength search, output model as model_crash.FILE_EXTENSION!");
        }

        if (config.get<IndexType>("useSeismogramTaper") == 2 || config.get<IndexType>("useSeismogramTaper") == 3) {                                                   
            seismogramTaper2D.apply(receivers.getSeismogramHandler()); 
        }
        seismogramTaper1D.apply(receivers.getSeismogramHandler());

        /* Normalize observed and synthetic data */
        if (config.get<bool>("normalizeTraces")) {
            receivers.getSeismogramHandler().normalize();
            receiversTrue.getSeismogramHandler().normalize();
        }

        misfitTest.setValue(shotIndTrue, dataMisfit.calc(receivers, receiversTrue, shotIndTrue));
    }

    commInterShot->sumArray(misfitTest.getLocalValues());

    HOST_PRINT(commAll, "\n======== Finished loop over test shots =========");
    HOST_PRINT(commAll, "\n================================================\n");

    return misfitTest.sum();
}

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
    if (configEM.getAndCatch("useRandomSource", 0) != 0) {  
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
    bool useStreamConfigEM = configEM.getAndCatch("useStreamConfig", false);
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
        if (configEM.getAndCatch("useRandomSource", 0) == 0) {  
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
            CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(configEM, sourceSettingsShot, *testmodelEM, modelCoordinatesEM, shotNumber);
        } else {
            CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(configEM, sourceSettingsShot, *testmodelPerShotEM, modelCoordinatesEM, shotNumber);
        }
        
        if (configEM.get<IndexType>("useReceiversPerShot") == 1) {
            receiversEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber);
            receiversTrueEM.init(configEM, modelCoordinatesEM, ctx, distEM, shotNumber);
        } else if (configEM.get<IndexType>("useReceiversPerShot") == 2) {
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
        if (configEM.get<IndexType>("useSeismogramTaper") == 2 || configEM.get<IndexType>("useSeismogramTaper") == 3) {
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

        if (configEM.get<IndexType>("useSeismogramTaper") == 2 || configEM.get<IndexType>("useSeismogramTaper") == 3) {                                                   
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

/*! \brief Initialize log-file
 *
 *
 \param comm Communicator
 \param logFilename Name of log-file
 */
template <typename ValueType>
void KITGPI::StepLengthSearch<ValueType>::initLogFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, std::string misfitType)
{
    int myRank = comm->getRank();
    if (myRank == MASTERGPI) {
        logFile.open(logFilename);
        logFile << "# Step length log file  \n";
        logFile << "# Misfit type = " << misfitType << "\n";
        logFile << "# Iteration 0 shows misfit of initial model of each workflow stage (only first two columns and last column is meaningful here)\n";
        logFile << "# Stage | Iteration | optimum step length | #Forward calculations | step length guess 1 | step length guess 2 | step length guess 3 | misfit of slg1 | misfit of slg2 | misfit of slg3 | final misfit of all shots\n";
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
void KITGPI::StepLengthSearch<ValueType>::appendToLogFile(scai::dmemo::CommunicatorPtr comm, IndexType workflowStage, scai::IndexType iteration, std::string logFilename, ValueType misfitSum)
{
    int myRank = comm->getRank();
    /* The following temporaries are only necessary because of a problem with LAMA: e.g. steplengthParabola.getValue(0).getValue<ValueType>() produces an error */
    ValueType steplengthParabola0 = steplengthParabola.getValue(0);
    ValueType steplengthParabola1 = steplengthParabola.getValue(1);
    ValueType steplengthParabola2 = steplengthParabola.getValue(2);
    ValueType misfitParabola0 = misfitParabola.getValue(0);
    ValueType misfitParabola1 = misfitParabola.getValue(1);
    ValueType misfitParabola2 = misfitParabola.getValue(2);

    if (myRank == MASTERGPI) {
        if (iteration == 0) {
            std::string filename(logFilename);
            logFile.open(filename, std::ios_base::app);
            logFile << std::scientific;
            logFile << std::setw(5) << workflowStage << std::setw(11) << iteration << std::setw(22) << 0 << std::setw(18) << 0 << std::setw(30) << 0 << std::setw(22) << 0 << std::setw(22) << 0 << std::setw(19) << 0 << std::setw(17) << 0 << std::setw(17) << 0 << std::setw(22) << misfitSum << "\n";
            logFile.close();
        } else {
            std::string filename(logFilename);
            logFile.open(filename, std::ios_base::app);
            logFile << std::scientific;
            logFile << std::setw(5) << workflowStage << std::setw(11) << iteration << std::setw(22) << steplengthOptimum << std::setw(18) << stepCalcCount << std::setw(30) << steplengthParabola0 << std::setw(22) << steplengthParabola1 << std::setw(22) << steplengthParabola2 << std::setw(19) << misfitParabola0 << std::setw(17) << misfitParabola1 << std::setw(17) << misfitParabola2 << std::setw(22) << misfitSum << "\n";
            logFile.close();
        }
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

template class KITGPI::StepLengthSearch<double>;
template class KITGPI::StepLengthSearch<float>;
