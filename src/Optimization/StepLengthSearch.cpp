#include "StepLengthSearch.hpp"
#include <string>
#include <iomanip>

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
void KITGPI::StepLengthSearch<ValueType>::run(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, ValueType steplengthInit, scai::lama::DenseVector<ValueType> currentMisfit)
{
    double start_t, end_t; /* For timing */
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    scai::dmemo::CommunicatorPtr comm = scai::dmemo::Communicator::getCommunicatorPtr();   // default communicator, set by environment variable SCAI_COMMUNICATOR
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                   // default context, set by environment variable SCAI_CONTEXT   
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    int testShotStart = config.get<int>("testShotStart");
    int testShotEnd = config.get<int>("testShotEnd");
    int testShotIncr = config.get<int>("testShotIncr");
    
    wavefields = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefields->init(ctx, dist);
    
    KITGPI::Acquisition::Seismogram<ValueType> synthetic(receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P));
    
    typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(KITGPI::Misfit::Factory<ValueType>::Create(config.get<std::string>("misfitType")));
    
    ValueType misfitTestSum;
    ValueType steplength;
    
    /* ------------------------------------------- */
    /* Set values for step length search           */
    /* ------------------------------------------- */
    stepCalcCount = 1;                                                           // number of forward calculations to find a proper (steplength, misfit) pair
    int maxStepCalc = config.get<int>("maxStepCalc");                            // maximum number of forward calculations to find a proper (steplength, misfit) pair 
    ValueType scalingFactor = config.get<ValueType>("scalingFactor");
    ValueType steplengthMin = config.get<ValueType>("steplengthMin");
    ValueType steplengthMax = config.get<ValueType>("steplengthMax");
    
    /* Save three pairs (steplength, misfit) for the parabolic fit */
    steplengthParabola.setSameValue(3, 0);
    misfitParabola.setSameValue(3, 0);

    /* ------------------------------------------------------------------------------------------------------ */
    /* Try to find a second step length which gives a smaller misfit and a third step length which gives again a larger misfit.
       The first step length is zero so the misfit from the gradient calculation can be used to save computation time.
       
       Separate 5 cases: 
       1) misfit 2 < misfit 1 AND misfit 3 > misfit 2: SL = minimum of parabola or SL = SLmax
       2) misfit 2 < misfit 1 AND misfit 3 < misfit 2: SL = SL3 or SL = SLmax 
       3) misfit 2 > misfit 1 AND misfit 3 < misfit 2 AND misfit 3 < misfit 1: SL = SL3 
       4) misfit 2 > misfit 1 AND misfit 3 < misfit 2 AND misfit 3 > misfit 1: SL = very small SL -> not used
       5) misfit 2 > misfit 1 AND misfit 3 > misfit 2: SL = very small SL */
    /* ------------------------------------------------------------------------------------------------------ */
    
    start_t = scai::common::Walltime::get();
    
    step2ok = false; // true if second step length decreases the misfit
    step3ok = false; // true if second step length decreases the misfit AND third step length increases the misfit relative to second step length */
    
    /* --- Save first step length (steplength=0 is used to save computational time) --- */
    misfitTestSum = 0;
    for(IndexType shotNumber = testShotStart ; shotNumber <= testShotEnd; shotNumber+=testShotIncr){
        misfitTestSum += currentMisfit.getValue(shotNumber);
    }
    misfitParabola.setValue(0, misfitTestSum);
   
    /* --- Save second step length (initial step length) in any case --- */
    HOST_PRINT(comm,"\nEstimation of 2nd steplength, forward test run no. " << stepCalcCount << " of maximum " << maxStepCalc <<"\n" );
    misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, receiversTrue, model, *wavefields, config, scaledGradient, *dataMisfit, steplengthInit);
    steplengthParabola.setValue(1, steplengthInit); 
    misfitParabola.setValue(1, misfitTestSum);  
    if ( misfitParabola.getValue(0) > misfitParabola.getValue(1) ){
        step2ok= true;}
    
    /* --- Search for a third step length - case: misfit was DECREASED --- */
    steplength = steplengthInit;
    while( step2ok == true && stepCalcCount < maxStepCalc ){
        HOST_PRINT(comm,"\nEstimation of 3rd steplength, forward test run no. " << stepCalcCount + 1 << " of maximum " << maxStepCalc <<"\n" );
        steplength *= scalingFactor;
        misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, receiversTrue, model, *wavefields, config, scaledGradient, *dataMisfit, steplength);
        steplengthParabola.setValue(2, steplength); 
        misfitParabola.setValue(2, misfitTestSum);
        if( misfitTestSum > misfitParabola.getValue(1)){
            step3ok= true;
            stepCalcCount += 1;
            break;}
        else{
            stepCalcCount += 1;}
    }
    
    /* --- Search for a third step length - case: misfit was INCREASED  --- */
    steplength = steplengthInit;
    while( step2ok == false && stepCalcCount < maxStepCalc ){
        HOST_PRINT(comm,"Estimation of 3rd steplength, forward test run no. " << stepCalcCount +1 << " of maximum " << maxStepCalc << "\n" );
        steplength /= scalingFactor;
        misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, receiversTrue, model, *wavefields, config, scaledGradient, *dataMisfit, steplength);
        steplengthParabola.setValue(2, steplength); 
        misfitParabola.setValue(2, misfitTestSum);
        if( misfitTestSum < misfitParabola.getValue(0)){
            step3ok= true;
            stepCalcCount += 1;
            break;}
        else{
            stepCalcCount += 1;}
    }
    
    /* Output of the three test step lengths and the corresponding misfit */
    for (int i = 0; i < 3; i++) {
        HOST_PRINT(comm,"\nSteplength " << i + 1 << ": "<<  steplengthParabola.getValue(i) << ", Corresponding misfit: " << misfitParabola.getValue(i) );
    }
    
    /* Set optimum step length */
    if( step2ok == true && step3ok == true ){
        steplengthOptimum = parabolicFit(steplengthParabola,misfitParabola);
        HOST_PRINT(comm,"\nApply parabolic fit");}
    else if( step2ok == true && step3ok == false ){
        steplengthOptimum = steplengthParabola.getValue(2);}
    else if( step2ok == false && step3ok == true ){
        steplengthOptimum = steplengthParabola.getValue(2);}
    else if( step2ok == false && step3ok == false ){
        steplengthOptimum = steplengthMin;
        HOST_PRINT(comm,"\nVariable steplengthMin used to update the model");}
        
    /* Check if accepted step length is smaller than maximally allowed step length */
    if ( steplengthOptimum > steplengthMax){
        steplengthOptimum = steplengthMax;
        HOST_PRINT(comm,"\nVariable steplengthMax used to update the model");}
        
    HOST_PRINT(comm,"\nOptimum step length: " << steplengthOptimum << "\n");    
        
    end_t = scai::common::Walltime::get();
    HOST_PRINT(comm,"\nFinished step length search in  " << end_t - start_t << " sec.\n\n\n" );

}

/*! \brief Parabolic fit
 *
 *
 \param xValues Vector of three values on the x-axis
 \param xValues Vector of three values on the y-axis
 */
template <typename ValueType>
ValueType KITGPI::StepLengthSearch<ValueType>::parabolicFit(scai::lama::DenseVector<ValueType> const &xValues,scai::lama::DenseVector<ValueType> const &yValues){
	
    SCAI_ASSERT(xValues.size() == 3, "xvalues must contain 3 values!");
    SCAI_ASSERT(yValues.size() == 3, "yvalues must contain 3 values!");    
    ValueType steplengthExtremum;
    scai::lama::DenseMatrix<ValueType> A(3,3);
    scai::lama::DenseMatrix<ValueType> invA;
    scai::lama::DenseVector<ValueType> coeff; // does the size need to be specified?
    scai::lama::DenseVector<ValueType> vectorOnes(3,1,scai::hmemo::Context::getContextPtr());
    scai::lama::DenseVector<ValueType> xValuesPow2 = xValues;
    xValuesPow2 *= xValues;
    
    A.setColumn(xValuesPow2, 0, scai::common::BinaryOp::COPY);
    A.setColumn(xValues, 1, scai::common::BinaryOp::COPY);
    A.setColumn(vectorOnes, 2, scai::common::BinaryOp::COPY);
    invA.invert(A); 
    
    coeff = invA * yValues;       
    steplengthExtremum = - coeff.getValue(1) / (2*coeff.getValue(0));
    
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
ValueType KITGPI::StepLengthSearch<ValueType>::calcMisfit(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, ValueType steplength)
{
    
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    scai::dmemo::DistributionPtr dist = wavefields.getRefVX().getDistributionPtr();
    scai::dmemo::CommunicatorPtr comm = scai::dmemo::Communicator::getCommunicatorPtr();   // default communicator, set by environment variable SCAI_COMMUNICATOR
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                   // default context, set by environment variable SCAI_CONTEXT   
        
//     double start_t, end_t; /* For timing */
    IndexType tStart = 0;
    IndexType tEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5); 
    std::string equationType = config.get<std::string>("equationType");
    int testShotStart = config.get<int>("testShotStart");
    int testShotEnd = config.get<int>("testShotEnd");
    int testShotIncr = config.get<int>("testShotIncr");
    
    scai::lama::DenseVector<ValueType> misfitTest(sources.getNumShots(), 0, ctx);
    
    // Implement a (virtual) copy constructor in the abstract base class to simplify the following code -> virtual constructor idiom!
    typename KITGPI::Modelparameter::Modelparameter<ValueType>::ModelparameterPtr testmodel(KITGPI::Modelparameter::Factory<ValueType>::Create(equationType));
    *testmodel = model;
    typename KITGPI::Gradient::Gradient<ValueType>::GradientPtr testgradient(KITGPI::Gradient::Factory<ValueType>::Create(equationType));
    *testgradient = scaledGradient;

    /* Update model */   
    *testgradient *= steplength;  
    *testmodel -= *testgradient;
        
    testmodel->prepareForModelling(config, ctx, dist, comm); 
    
    // later it should be possible to select only a subset of shots for the step length search
    for (IndexType shotNumber = testShotStart ; shotNumber <= testShotEnd; shotNumber+=testShotIncr) {
        
        wavefields.resetWavefields();

        HOST_PRINT(comm, "\n=============== Shot " << shotNumber + 1 << " of " << sources.getNumShots() << " ===================\n");
        HOST_PRINT(comm, "\n--------------- Start Test Forward -------------------\n");

        sources.init(config, ctx, dist, shotNumber);
        
//         start_t = scai::common::Walltime::get();
        solver.run(receivers, sources, *testmodel, wavefields, derivatives, tStart, tEnd, config.get<ValueType>("DT"));
//         end_t = scai::common::Walltime::get();
//         HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");
        
        receiversTrue.getSeismogramHandler().readFromFileRaw(config.get<std::string>("FieldSeisName")  + ".shot_" + std::to_string(shotNumber) + ".mtx", 1);
        
        misfitTest.setValue(shotNumber, dataMisfit.calc(receivers, receiversTrue));
        
    }
    
    HOST_PRINT(comm, "\n======== Finished loop over test shots =========");
    HOST_PRINT(comm, "\n================================================\n");
    
    return misfitTest.sum();
            
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
        if(iteration==0){
            std::string filename(logFilename);
            logFile.open(filename, std::ios_base::app);
            logFile <<  std::scientific ;
            logFile << std::setw(5) << workflowStage << std::setw(11) << iteration << std::setw(22) << 0 << std::setw(18) << 0 << std::setw(30) << 0 << std::setw(22) << 0 << std::setw(22) << 0 << std::setw(19) << 0 << std::setw(17) << 0 << std::setw(17) << 0 << std::setw(22) << misfitSum << "\n" ;
            logFile.close();
        } else {
            std::string filename(logFilename);
            logFile.open(filename, std::ios_base::app);
            logFile <<  std::scientific ;
            logFile << std::setw(5) << workflowStage << std::setw(11) << iteration << std::setw(22) << steplengthOptimum << std::setw(18) << stepCalcCount << std::setw(30) << steplengthParabola0 << std::setw(22) << steplengthParabola1 << std::setw(22) << steplengthParabola2 << std::setw(19) << misfitParabola0 << std::setw(17) << misfitParabola1 << std::setw(17) << misfitParabola2 << std::setw(22) << misfitSum << "\n" ;
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

