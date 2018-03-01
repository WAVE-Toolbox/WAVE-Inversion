#include "StepLengthSearch.hpp"
#include <string>
#include <iomanip>

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
 \param steplength_init Initial steplength
 \param currentMisfit Current misfit
 */
template <typename ValueType>
void StepLengthSearch<ValueType>::run(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, scai::lama::Scalar steplength_init, scai::lama::Scalar currentMisfit)
{
    double start_t, end_t; /* For timing */
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    scai::dmemo::CommunicatorPtr comm = scai::dmemo::Communicator::getCommunicatorPtr();   // default communicator, set by environment variable SCAI_COMMUNICATOR
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                   // default context, set by environment variable SCAI_CONTEXT   
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    
    wavefields = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefields->init(ctx, dist);
    
    KITGPI::Acquisition::Seismogram<ValueType> synthetic(receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P));
    
    typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr dataMisfit(KITGPI::Misfit::Factory<ValueType>::Create(config.get<std::string>("misfitType")));
    
    scai::lama::Scalar misfitTestSum;
    scai::lama::Scalar steplength;
    
    /* ------------------------------------------- */
    /* Set values for step length search           */
    /* ------------------------------------------- */
    stepCalcCount = 0;                                                       // number of calculations to find a proper (steplength, misfit) pair
    int maxStepCalc = config.get<int>("MaxStepCalc");                            // maximum number of calculations to find a proper (steplength, misfit) pair 
    scai::lama::Scalar scalingFactor = config.get<ValueType>("scalingFactor");
    scai::lama::Scalar steplengthMin = config.get<ValueType>("steplengthMin");
    scai::lama::Scalar steplengthMax = config.get<ValueType>("steplengthMax");
    
    /* Save three pairs (steplength, misfit) for the parabolic fit */
    steplengthParabola.allocate(3);
    misfitParabola.allocate(3);
    steplengthParabola.assign(0.0);
    misfitParabola.assign(0.0);
    misfitParabola.setValue(0, currentMisfit);  // use pair (0, current misfit) to save computation time

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
   
    /* --- Save second step length in any case --- */
    HOST_PRINT(comm,"\nEstimation of 2nd steplength \n\n" );
    misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, receiversTrue, model, *wavefields, config, scaledGradient, *dataMisfit, steplength_init);
    steplengthParabola.setValue(1, steplength_init); 
    misfitParabola.setValue(1, misfitTestSum);  
    if ( misfitParabola.getValue(0) > misfitParabola.getValue(1) ){
        step2ok= true;}
    
    /* --- Search for a third step length - case: misfit was DECREASED --- */
    steplength = steplength_init;
    while( step2ok == true && stepCalcCount < maxStepCalc ){
	HOST_PRINT(comm,"Estimation of 3rd steplength no. " << stepCalcCount +1 << " of " << maxStepCalc <<"\n" );
        steplength *= scalingFactor;
        misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, receiversTrue, model, *wavefields, config, scaledGradient, *dataMisfit, steplength);
        steplengthParabola.setValue(2, steplength); 
        misfitParabola.setValue(2, misfitTestSum);
        if( misfitTestSum > misfitParabola.getValue(1)){
            step3ok= true;
            break;}
        else{
            stepCalcCount += 1;}
    }
    
    /* --- Search for a third step length - case: misfit was INCREASED  --- */
    steplength = steplength_init;
    while( step2ok == false && stepCalcCount < maxStepCalc ){
	HOST_PRINT(comm,"Estimation of 3rd steplength no. " << stepCalcCount+1 << " of " << maxStepCalc <<"\n" );
        steplength /= scalingFactor;
        misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, receiversTrue, model, *wavefields, config, scaledGradient, *dataMisfit, steplength);
        steplengthParabola.setValue(2, steplength); 
        misfitParabola.setValue(2, misfitTestSum);
        if( misfitTestSum < misfitParabola.getValue(0)){
            step3ok= true;
            break;}
        else{
            stepCalcCount += 1;}
    }
    
    /* Set optimum step length */
    if( step2ok == true && step3ok == true ){
        HOST_PRINT(comm,"Apply parabolic fit\n" );
        steplengthOptimum = parabolicFit(steplengthParabola,misfitParabola);}
    else if( step2ok == true && step3ok == false ){
        steplengthOptimum = steplengthParabola.getValue(2);}
    else if( step2ok == false && step3ok == true ){
        steplengthOptimum = steplengthParabola.getValue(2);}
    else if( step2ok == false && step3ok == false ){
        steplengthOptimum = steplengthMin;}
        
        
    /* Check if accepted step length is smaller than maximally allowed step length */
    if ( steplengthOptimum > steplengthMax){
        steplengthOptimum = steplengthMax;}
        
    end_t = scai::common::Walltime::get();
    HOST_PRINT(comm,"\nFinished step length search in  " << end_t - start_t << " sec.\n\n\n" );
//     for (int i = 0; i < 3; i++) {
//         HOST_PRINT(comm,"Steplength " << i << ": "<<  steplengthParabola.getValue(i) << ", Corresponding misfit: " << misfitParabola.getValue(i) << std::endl);
//     }
//     HOST_PRINT(comm,"Accepted step length: " << steplengthOptimum << ", Corresponding misfit: " << misfitTestSum << std::endl);
    
}

/*! \brief Parabolic fit
 *
 *
 \param xValues Vector of three values on the x-axis
 \param xValues Vector of three values on the y-axis
 */
template <typename ValueType>
scai::lama::Scalar StepLengthSearch<ValueType>::parabolicFit(scai::lama::DenseVector<ValueType> const &xValues,scai::lama::DenseVector<ValueType> const &yValues){
	
    SCAI_ASSERT(xValues.size() == 3, "xvalues must contain 3 values!");
    SCAI_ASSERT(yValues.size() == 3, "yvalues must contain 3 values!");
    scai::lama::Scalar steplengthExtremum;
    scai::lama::DenseMatrix<ValueType> A(3,3);
    scai::lama::DenseMatrix<ValueType> invA;
    scai::lama::DenseVector<ValueType> coeff; // does the size need to be specified?
    scai::lama::DenseVector<ValueType> vectorOnes(3,1,scai::hmemo::Context::getContextPtr());
    scai::lama::DenseVector<ValueType> xValuesPow2 = xValues;
    xValuesPow2 *= xValues;
    
    A.setColumn(xValuesPow2, 0, scai::common::binary::BinaryOp::COPY);
    A.setColumn(xValues, 1, scai::common::binary::BinaryOp::COPY);
    A.setColumn(vectorOnes, 2, scai::common::binary::BinaryOp::COPY);
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
scai::lama::Scalar StepLengthSearch<ValueType>::calcMisfit(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, scai::lama::Scalar steplength)
{
    
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    scai::dmemo::DistributionPtr dist = wavefields.getRefVX().getDistributionPtr();
    scai::dmemo::CommunicatorPtr comm = scai::dmemo::Communicator::getCommunicatorPtr();   // default communicator, set by environment variable SCAI_COMMUNICATOR
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                   // default context, set by environment variable SCAI_CONTEXT   
        
   // double start_t, end_t; /* For timing */
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
        
        wavefields.reset();
        sources.init(config, ctx, dist, shotNumber);
        
     //   start_t = scai::common::Walltime::get();
        solver.run(receivers, sources, *testmodel, wavefields, derivatives, tStart, tEnd, config.get<ValueType>("DT"));
       // end_t = scai::common::Walltime::get();
        //HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");
        
        receiversTrue.getSeismogramHandler().readFromFileRaw(config.get<std::string>("FieldSeisName")  + ".shot_" + std::to_string(shotNumber) + ".mtx", 1);
        
        misfitTest.setValue(shotNumber, dataMisfit.calc(receivers, receiversTrue));
        
    }
    
    return misfitTest.sum();
            
}

/*! \brief Initialize log-file
 *
 *
 \param comm Communicator
 \param logFilename Name of log-file
 */
template <typename ValueType>
void StepLengthSearch<ValueType>::initLogFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, std::string misfitType)
{
    int myRank = comm->getRank(); 
    if (myRank == MASTERGPI) {
        logFile.open(logFilename);
        logFile << "# Step length log file  \n";
        logFile << "# Misfit type = " << misfitType << "\n";
        logFile << "# Iteration 0 shows misfit of initial model (only first and last column is meaningful here)\n";
        logFile << "# Iteration | optimum step length | #Forward calculations | step length guess 1 | step length guess 2 | step length guess 3 | misfit of slg1 | misfit of slg2 | misfit of slg3 | final misfit of all shots\n";
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
void StepLengthSearch<ValueType>::appendToLogFile(scai::dmemo::CommunicatorPtr comm, IndexType iteration, std::string logFilename, scai::lama::Scalar misfitSum)
{
    int myRank = comm->getRank(); 
    /* The following temporaries are only necessary because of a problem with LAMA: e.g. steplengthParabola.getValue(0).getValue<ValueType>() produces an error */
    scai::lama::Scalar steplengthParabola0 = steplengthParabola.getValue(0);
    scai::lama::Scalar steplengthParabola1 = steplengthParabola.getValue(1);
    scai::lama::Scalar steplengthParabola2 = steplengthParabola.getValue(2);
    scai::lama::Scalar misfitParabola0 = misfitParabola.getValue(0);
    scai::lama::Scalar misfitParabola1 = misfitParabola.getValue(1);
    scai::lama::Scalar misfitParabola2 = misfitParabola.getValue(2);
    
    if (myRank == MASTERGPI) {
        std::string filename(logFilename);
        logFile.open(filename, std::ios_base::app);
        logFile <<  std::scientific ;
        logFile << std::setw(8) << iteration << std::setw(22) << steplengthOptimum.getValue<ValueType>() << std::setw(18) << stepCalcCount << std::setw(30) << steplengthParabola0.getValue<ValueType>() << std::setw(22) << steplengthParabola1.getValue<ValueType>() << std::setw(22) << steplengthParabola2.getValue<ValueType>() << std::setw(19) << misfitParabola0.getValue<ValueType>() << std::setw(17) << misfitParabola1.getValue<ValueType>() << std::setw(17) << misfitParabola2.getValue<ValueType>() << std::setw(22) << misfitSum.getValue<ValueType>() << "\n" ;
        logFile.close();
    }                             
}

/*! \brief Get optimum steplength
 *
 *
 */
template <typename ValueType>
scai::lama::Scalar const &StepLengthSearch<ValueType>::getSteplength()
{
    return (steplengthOptimum);
}


template class StepLengthSearch<double>;
template class StepLengthSearch<float>;

