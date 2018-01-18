#include "StepLengthSearch.hpp"

template <typename ValueType>
void StepLengthSearch<ValueType>::calc(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, scai::lama::Scalar steplength_init, scai::lama::Scalar currentMisfit)
{
    
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    scai::dmemo::CommunicatorPtr comm = scai::dmemo::Communicator::getCommunicatorPtr();   // default communicator, set by environment variable SCAI_COMMUNICATOR
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                   // default context, set by environment variable SCAI_CONTEXT   
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    
    wavefields = KITGPI::Wavefields::Factory<ValueType>::Create(dimension, equationType);
    wavefields->init(ctx, dist);
    
    KITGPI::Acquisition::Seismogram<ValueType> truedata(receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P));
    KITGPI::Acquisition::Seismogram<ValueType> synthetic(receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P));
    
    scai::lama::Scalar misfitTestSum;
    scai::lama::Scalar steplength;
    
    int maxStepCalc = 4;      // maximum number of calculations to find a proper (steplength, misfit) pair 
    int stepCalcCount = 0;    // number of calculations to find a proper (steplength, misfit) pair
    scai::lama::Scalar scalingFactor = 2;
    scai::lama::Scalar steplengthMin = 0.001; 
    scai::lama::Scalar steplengthMax = 0.1; // 0.1 shows minimum msifit of ~97
        
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
    
    HOST_PRINT(comm,"Start step length search\n\n" );
    
    step2ok = false; // true if second step length decreases the misfit
    step3ok = false; // true if second step length decreases the misfit AND third step length increases the misfit relative to second step length */
   
    /* --- Save first step length in any case --- */
    misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, model, *wavefields, config, scaledGradient, steplength_init);
    steplengthParabola.setValue(1, steplength_init); 
    misfitParabola.setValue(1, misfitTestSum);  
    if ( misfitParabola.getValue(0) > misfitParabola.getValue(1) ){
        step2ok= true;}
    
    /* --- Search for a second step length - case: misfit was DECREASED --- */
    steplength = steplength_init;
    while( step2ok == true && stepCalcCount < maxStepCalc ){
        steplength *= scalingFactor;
        misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, model, *wavefields, config, scaledGradient, steplength);
        steplengthParabola.setValue(2, steplength); 
        misfitParabola.setValue(2, misfitTestSum);
        if( misfitTestSum > misfitParabola.getValue(1)){
            step3ok= true;
            break;}
        else{
            stepCalcCount += 1;}
    }
    
    /* --- Search for a second step length - case: misfit was INCREASED  --- */
    steplength = steplength_init;
    while( step2ok == false && stepCalcCount < maxStepCalc ){
        steplength /= scalingFactor;
        misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, model, *wavefields, config, scaledGradient, steplength);
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
        HOST_PRINT(comm,"Apply parabolic fit\n\n" );
        this->parabolicFit();
        steplengthOptimum = steplengthExtremum;}
    else if( step2ok == true && step3ok == false ){
        steplengthOptimum = steplengthParabola.getValue(2);}
    else if( step2ok == false && step3ok == true ){
        steplengthOptimum = steplengthParabola.getValue(2);}
    else if( step2ok == false && step3ok == false ){
        steplengthOptimum = steplengthMin;}
        
        
    /* Check if accepted step length is smaller than maximally allowed step length */
    if ( steplengthOptimum > steplengthMax){
        steplengthOptimum = steplengthMax;}
    
    HOST_PRINT(comm,"Finished step length search\n\n" );
    for (int i = 0; i < 3; i++) {
        HOST_PRINT(comm,"Steplength " << i << ": "<<  steplengthParabola.getValue(i) << ", Corresponding misfit: " << misfitParabola.getValue(i) << std::endl);
    }
    HOST_PRINT(comm,"Accepted step length: " << steplengthOptimum << ", Corresponding misfit: " << misfitTestSum << std::endl);
    
    /* -------------------------------------- */
    /* Find second step length (out of three) */
    /* -------------------------------------- */
//     while( step2ok == false ){
//         
//         if(stepCalcCount >= maxStepCalc)
//             break;
//         
//         misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, model, *wavefields, config, scaledGradient, steplength);
//         misfitParabola.setValue(1, misfitTestSum);
//         steplengthParabola.setValue(1, steplength);
//         if(misfitTestSum < currentMisfit){
//             step2ok = true; // use set method?
//             steplength = scalingFactor * steplength;}
//         else {
//             steplength /= scalingFactor;
//             stepCalcCount += 1;}
//     }
    
    /* ----------------------------------------------------------------------------------- */
    /* Find third step length in case that the second step length DOES decrease the misfit */
    /* ----------------------------------------------------------------------------------- */
//     stepCalcCount = 0;
//     
//     if( step2ok == true ){
//         
//         while( step3ok == false ){
//             if(stepCalcCount >= maxStepCalc){
//                 if( steplength > steplengthMax){
//                     steplengthOptimum = steplengthMax;
//                     return;}
//                 break;} 
//         
//             misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, model, *wavefields, config, scaledGradient, steplength);
//             misfitParabola.setValue(2, misfitTestSum);
//             steplengthParabola.setValue(2, steplength);
//             if(misfitTestSum > misfitParabola.getValue(1)){
//                 step3ok = true; // use set method?
//                 stepCalcCount = 0;}  // do not set stepCalcCount to zero??
//             else if( misfitTestSum < misfitParabola.getValue(1) ){
//                 steplength = scalingFactor * steplength;
//                 stepCalcCount += 1;}
//         }
//     }
//     
//     /* Fit parabola */
//     if( step2ok == true && step3ok == true ){
//         this->parabolicFit();
//         steplengthOptimum = steplengthExtremum;
//         return;
//     }
    
    /* --------------------------------------------------------------------------------------- */
    /* Find third step length in case that the second step length DOES NOT decrease the misfit */
    /* --------------------------------------------------------------------------------------- */
//     stepCalcCount = 0;
//     steplength = steplength_init;
//     
//     if( step2ok == false ){
//         
//         while(stepCalcCount < maxStepCalc){
//             
//             misfitTestSum = this->calcMisfit(solver, derivatives, receivers, sources, model, *wavefields, config, scaledGradient, steplength);
//             misfitParabola.setValue(2, misfitTestSum);
//             steplengthParabola.setValue(2, steplength);
//             if(misfitTestSum < misfitParabola.getValue(0)){
//                 stepCalcCount = 0;
//                 steplengthOptimum = steplength; // set condition with steplengthMax?
//                 return;}
//             else{
//                 steplength = scalingFactor * steplength;
//                 }
//         }
//         
//         steplengthOptimum = steplengthMin;
//         return;
//     }
    
    
}


template <typename ValueType>
void StepLengthSearch<ValueType>::parabolicFit(){
      
    scai::lama::DenseMatrix<ValueType> A(3,3);
    scai::lama::DenseMatrix<ValueType> invA;
    scai::lama::DenseVector<ValueType> coeffParabola; // does the size need to be specified?
    scai::lama::DenseVector<ValueType> vectorOnes(3,1,scai::hmemo::Context::getContextPtr());
    scai::lama::DenseVector<ValueType> steplengthParabolaPow2 = steplengthParabola;
    steplengthParabolaPow2 *= steplengthParabola;
    
    A.setColumn(steplengthParabolaPow2, 0, scai::common::binary::BinaryOp::COPY);
    A.setColumn(steplengthParabola, 1, scai::common::binary::BinaryOp::COPY);
    A.setColumn(vectorOnes, 2, scai::common::binary::BinaryOp::COPY);
    invA.invert(A); 
    
    coeffParabola = invA * misfitParabola;       
    steplengthExtremum = - coeffParabola.getValue(1) / (2*coeffParabola.getValue(0));
    
}


template <typename ValueType>
scai::lama::Scalar StepLengthSearch<ValueType>::calcMisfit(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, KITGPI::Gradient::Gradient<ValueType> &scaledGradient, scai::lama::Scalar steplength)
{
    
    /* ------------------------------------------- */
    /* Get distribution, communication and context */
    /* ------------------------------------------- */
    scai::dmemo::DistributionPtr dist = wavefields.getRefVX().getDistributionPtr();
    scai::dmemo::CommunicatorPtr comm = scai::dmemo::Communicator::getCommunicatorPtr();   // default communicator, set by environment variable SCAI_COMMUNICATOR
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                   // default context, set by environment variable SCAI_CONTEXT   
        
    double start_t, end_t; /* For timing */
    IndexType tStart = 0;
    IndexType tEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5); 
    std::string equationType = config.get<std::string>("equationType");

    KITGPI::Acquisition::Seismogram<ValueType> truedata(receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P));
    KITGPI::Acquisition::Seismogram<ValueType> synthetic(receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P));
    
    scai::lama::Scalar misfitTestShot;
    scai::lama::Scalar misfitTestSum;
    
    // Implement a (virtual) copy constructor in the abstract base class to simplify/generalize the following code -> virtual constructor idiom!
    typename KITGPI::Modelparameter::Modelparameter<ValueType>::ModelparameterPtr testmodel(KITGPI::Modelparameter::Factory<ValueType>::Create(equationType));
    testmodel->init(config, ctx, dist);
    testmodel->setVelocityP(model.getVelocityP());
    typename KITGPI::Gradient::Gradient<ValueType>::GradientPtr testgradient(KITGPI::Gradient::Factory<ValueType>::Create(equationType));
    testgradient->init(ctx, dist);
    testgradient->setVelocityP(scaledGradient.getVelocityP());

    /* Update model */   
    *testgradient *= steplength;
    *testmodel -= *testgradient;
        
    testmodel->prepareForModelling(config, ctx, dist, comm); 
    
    misfitTestSum = 0;
    
    std::string fieldSeisName("seismograms/rectangle.true");
    
    // later it should be possible to select only a subset of shots for the step length search
    for (IndexType shotNumber = 0; shotNumber < sources.getNumShots(); shotNumber++) {
        
        misfitTestShot = 0;
        wavefields.reset();
        sources.init(config, ctx, dist, shotNumber);
        
        start_t = scai::common::Walltime::get();
        solver.run(receivers, sources, *testmodel, wavefields, derivatives, tStart, tEnd, config.get<ValueType>("DT"));
        end_t = scai::common::Walltime::get();
        HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");
        
        truedata.readFromFileRaw(fieldSeisName +".It0" + ".shot" + std::to_string(shotNumber) + ".p.mtx", receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P).getData().getRowDistributionPtr(), NULL);
        synthetic = receivers.getSeismogramHandler().getSeismogram(KITGPI::Acquisition::SeismogramType::P);
        
        synthetic -= truedata;
        misfitTestShot = 0.5*synthetic.getData().l2Norm();   // misfit of one shot
        misfitTestSum += misfitTestShot;                     // misfit sum of selected shots
        
    }
    
    return misfitTestSum;
            
}

template <typename ValueType>
scai::lama::Scalar const &StepLengthSearch<ValueType>::getSteplength()
{
    return (steplengthOptimum);
}


template class StepLengthSearch<double>;
template class StepLengthSearch<float>;

