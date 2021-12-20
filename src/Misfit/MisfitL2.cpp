#include "MisfitL2.hpp"

/*! \brief init function
 *
 \param config Configuration handle 
 \param misfitTypeHistory a history vector to record the used times for each misfitType 
 \param numshots number of shots 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::init(KITGPI::Configuration::Configuration config, std::vector<scai::IndexType> misfitTypeHistory, scai::IndexType numshots)
{    
    // transform to lower cases
    misfitType = config.get<std::string>("misfitType");
    multiMisfitType = config.getAndCatch("multiMisfitType", multiMisfitType);
    std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::tolower);
    saveMultiMisfits = config.getAndCatch("saveMultiMisfits", false);
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                 // default context, set by environment variable SCAI_CONTEXT 
    scai::lama::DenseVector<ValueType> temp1(numshots, 0, ctx);  
    misfitTypeShots = temp1; // initialize misfitTypeShots to 0 because sumArray will be applied on it when misfitType=l2567892
    numMisfitTypes = multiMisfitType.length()-1;
    uniqueMisfitTypes.clear();
    for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) {  
        std::string temp = multiMisfitType.substr(iMisfitType+1, 1);
        uniqueMisfitTypes.push_back(atoi(temp.c_str()));
    }
    if (saveMultiMisfits || misfitType.length() > 2) {
        scai::lama::DenseVector<ValueType> misfitPerIt(numshots, 0, ctx);
        for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) {  
            misfitStorageL2.push_back(misfitPerIt); 
        }
        scai::IndexType misfitStorageL2Size = misfitStorageL2.size();
        if (misfitStorageL2Size == numMisfitTypes) {
            misfitSum0Ratio.clear();
            scai::lama::DenseVector<ValueType> temp3(numshots, 1, ctx);
            for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) {  
                misfitSum0Ratio.push_back(temp3); 
            }
        }
    } 
    if (misfitType.length() == 2) {
        std::string temp = misfitType.substr(1, 1);
        scai::IndexType misfitTypeNo = atoi(temp.c_str());
        misfitTypeShots = misfitTypeNo;
        if (!saveMultiMisfits) {
            misfitSum0Ratio.clear();
            scai::lama::DenseVector<ValueType> temp3(numshots, 1, ctx);
            for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) {  
                misfitSum0Ratio.push_back(temp3); 
            }
        }
    }
    if (misfitType.length() > 2 && misfitType.substr(numMisfitTypes+1, 1).compare("1") == 0) {
        scai::IndexType randMisfitTypeInd;
        scai::IndexType maxiterations = config.get<scai::IndexType>("maxIterations");
        std::srand((int)time(0));
        scai::IndexType maxcount = maxiterations * numshots / uniqueMisfitTypes.size() * 1.2;
        for (int shotInd = 0; shotInd < numshots; shotInd++) {         
            randMisfitTypeInd = std::rand() % uniqueMisfitTypes.size();
            if (misfitTypeHistory.at(randMisfitTypeInd) >= maxcount) {
                shotInd--;
            } else {
                misfitTypeShots.setValue(shotInd, uniqueMisfitTypes.at(randMisfitTypeInd));
                misfitTypeHistory.at(randMisfitTypeInd)++;
            }
        }
    }
}

/*! \brief Append to misfitType-file
*
\param comm Communicator
\param logFilename Name of log-file
\param stage inversion stage
\param iteration inversion iteration
*/
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::appendMisfitTypeShotsToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration)
{      
    int myRank = comm->getRank();  
    if (saveMultiMisfits && myRank == MASTERGPI) {
        std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::toupper);
        std::ofstream outputFile; 
        std::string misfitTypeFilename = logFilename.substr(0, logFilename.length()-4) + ".misfitType" + logFilename.substr(logFilename.length()-4, 4);
        if (stage == 1 && iteration == 0) {
            outputFile.open(misfitTypeFilename);
            outputFile << "# Misfit type records during inversion\n"; 
            outputFile << "# Misfit type = " << misfitType << "\n"; 
            outputFile << "# Stage | Iteration | misfitTypes\n"; 
        } else {                    
            outputFile.open(misfitTypeFilename, std::ios_base::app);
            outputFile << std::scientific;
        }
        outputFile << std::setw(5) << stage << std::setw(10) << iteration;
        for (int shotInd = 0; shotInd < misfitTypeShots.size(); shotInd++) { 
            if (shotInd == 0) {
                outputFile << std::setw(9) << (int) misfitTypeShots.getValue(shotInd);
            } else {
                outputFile << std::setw(4) << (int) misfitTypeShots.getValue(shotInd);
            }
        }
        outputFile << "\n";
        outputFile.close();
        std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::tolower);
    }
}

/*! \brief Append to misfitType-file
*
\param comm Communicator
\param logFilename Name of log-file
\param stage inversion stage
\param iteration inversion iteration
*/
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::appendMisfitPerShotToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration)
{      
    int myRank = comm->getRank();  
    if (saveMultiMisfits && myRank == MASTERGPI) {
        std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::toupper);
        std::ofstream outputFile; 
        std::string misfitTypeFilename = logFilename.substr(0, logFilename.length()-4) + ".pershot" + logFilename.substr(logFilename.length()-4, 4);
        if (stage == 1 && iteration == 0) {
            outputFile.open(misfitTypeFilename);
            outputFile << "# Misfit per shot during inversion\n"; 
            outputFile << "# Misfit type = " << misfitType << "\n"; 
            outputFile << "# Stage | Iteration | misfits\n"; 
        } else {                    
            outputFile.open(misfitTypeFilename, std::ios_base::app);
            outputFile << std::scientific;
        }
        outputFile << std::setw(5) << stage << std::setw(10) << iteration;
        scai::IndexType misfitStorageL2Size = misfitStorageL2.size();
        scai::lama::DenseVector<ValueType> misfitPerIt = misfitStorageL2.at(misfitStorageL2Size - numMisfitTypes);
        for (int shotInd = 0; shotInd < misfitTypeShots.size(); shotInd++) { 
            if (shotInd == 0) {
                outputFile << std::setw(9) << misfitPerIt.getValue(shotInd);
            } else {
                outputFile << std::setw(4) << misfitPerIt.getValue(shotInd);
            }
        }
        outputFile << "\n";
        outputFile.close();
        std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::tolower);
    }
}

/*! \brief Append to misfits-file
*
\param comm Communicator
\param logFilename Name of log-file
\param stage inversion stage
\param iteration inversion iteration
*/
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::appendMisfitsToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration)
{      
    int myRank = comm->getRank();  
    if (saveMultiMisfits && myRank == MASTERGPI) {
        std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::toupper);
        std::ofstream outputFile; 
        std::string misfitTypeFilename = logFilename.substr(0, logFilename.length()-3) + multiMisfitType.substr(0, numMisfitTypes+1) + logFilename.substr(logFilename.length()-4, 4);
        if (stage == 1 && iteration == 0) {
            outputFile.open(misfitTypeFilename);
            outputFile << "# Misfit per type during inversion\n"; 
            outputFile << "# Misfit type = " << misfitType << "\n"; 
            outputFile << "# Stage | Iteration"; 
            for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) { 
                outputFile << " | " << std::setw(7) << uniqueMisfitTypes.at(iMisfitType) << std::setw(7);
            }
            outputFile << "\n"; 
        } else {                    
            outputFile.open(misfitTypeFilename, std::ios_base::app);
            outputFile << std::scientific;
        }
        outputFile << std::setw(5) << stage << std::setw(10) << iteration;
        
        scai::IndexType misfitStorageL2Size = misfitStorageL2.size();
        for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) { 
            scai::lama::DenseVector<ValueType> misfitPerIt = misfitStorageL2.at(misfitStorageL2Size - numMisfitTypes + iMisfitType);
            if (iMisfitType == 0) {
                outputFile << std::setw(18) << misfitPerIt.sum();
            } else {
                outputFile << std::setw(14) << misfitPerIt.sum();
            }
        }
        outputFile << "\n";
        outputFile.close(); 
        std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::tolower);
    }
}

/*! \brief Append to misfits-file
*
\param comm Communicator
\param logFilename Name of log-file
\param stage inversion stage
\param iteration inversion iteration
*/
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot)
{      
    if (saveMultiMisfits || misfitType.length() > 2) {
        for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) { 
            scai::IndexType misfitStorageL2Size = misfitStorageL2.size();
            commInterShot->sumArray(misfitStorageL2.at(misfitStorageL2Size - numMisfitTypes + iMisfitType).getLocalValues());   
        }
        if (misfitType.length() > 2 && misfitType.substr(numMisfitTypes+1, 1).compare("2") == 0) {
            commInterShot->sumArray(misfitTypeShots.getLocalValues()); 
        }
//         std::cout<< "misfitTypeShots sumArray : ";
//         for (int i = 0; i < misfitTypeShots.size(); i++) { 
//             std::cout<< misfitTypeShots.getValue(i) << " ";
//         }
//         std::cout<< std::endl;
    }
}
/*! \brief Return the L2-norm of seismograms stored in the given receiver objects (note: the misfit is summed over all components, i.e. vx,vy,vz,p)
 *
 *
 \param receiversSyn Receiver object which stores the synthetic data 
 \param receiversObs Receiver object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs, scai::IndexType shotInd)
{        
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObs;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerSyn;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerObs;
    
    ValueType misfit = 0;
    ValueType misfitSum = 0;
    scai::IndexType count = 0;

    seismoHandlerSyn = receiversSyn.getSeismogramHandler();
    seismoHandlerObs = receiversObs.getSeismogramHandler();
          
    /* Note that the misfit of different components (p,vx,vy,vz) is summed up. If p and v? is used at the same time, this could cause problems because they could have different scales.
       For different velocity components it should be ok. */
    
    scai::IndexType misfitStorageL2Size = 0;
    scai::IndexType misfitTypeShotL2 = misfitTypeShots.getValue(shotInd);
    bool errorMisfitType = false;
    numMisfitTypes = misfitSum0Ratio.size();
    for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) {  
        if (misfitTypeShotL2 == 0 || misfitTypeShotL2 == uniqueMisfitTypes.at(iMisfitType)) {
            errorMisfitType = true;
            break;
        }        
    }
    SCAI_ASSERT_ERROR(errorMisfitType, "misfitType = L" + std::to_string(misfitTypeShotL2));
        
    for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) {  
        if (saveMultiMisfits || misfitType.length() > 2)
            misfitTypeShotL2 = uniqueMisfitTypes.at(iMisfitType);
            
        for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
            seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
            seismogramObs = seismoHandlerObs.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
            switch (misfitTypeShotL2) {
            case 5:
                misfit = this->calcL2EnvelopeWeighted(seismogramSyn, seismogramObs);
                break;
            case 6:
                misfit = this->calcL2AGC(seismogramSyn, seismogramObs);
                break;
            case 7:
                misfit = this->calcL2Normalized(seismogramSyn, seismogramObs);
                break;
            case 8:
                misfit = this->calcL2Envelope(seismogramSyn, seismogramObs);
                break;
            case 9:
                misfit = this->calcL2InstantaneousPhase(seismogramSyn, seismogramObs);
                break;
            default:
                misfit = this->calcL2(seismogramSyn, seismogramObs);
                break;
            } 
            misfitSum += misfit;
            if (misfit != 0) count++;
        }  
        
        if (saveMultiMisfits || misfitType.length() > 2) {
            misfitStorageL2Size = misfitStorageL2.size();
            if (misfitStorageL2Size == numMisfitTypes && iMisfitType > 0) {
                misfitSum0Ratio.at(iMisfitType).setValue(shotInd, misfitStorageL2.at(0).getValue(shotInd) / (misfitSum/count/misfitTypeShots.size()));
            }
            misfitStorageL2.at(misfitStorageL2Size - numMisfitTypes + iMisfitType).setValue(shotInd, misfitSum/count/misfitTypeShots.size()*misfitSum0Ratio.at(iMisfitType).getValue(shotInd));
            
            if (iMisfitType < numMisfitTypes - 1) {
                misfitSum = 0;
                count = 0;
            }
        } else {      
            break;
        }
    }
    if (misfitType.length() > 2 && misfitType.substr(numMisfitTypes+1, 1).compare("2") == 0) {
        scai::IndexType misfitRatioMaxInd = 0;
        if (misfitStorageL2Size > numMisfitTypes) { // iteration > 0
            std::vector<ValueType> misfitRatio(numMisfitTypes, 0);  
            for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) {           
                misfitRatio.at(iMisfitType) = (misfitStorageL2.at(misfitStorageL2Size - 2 * numMisfitTypes + iMisfitType).getValue(shotInd) - misfitStorageL2.at(misfitStorageL2Size - numMisfitTypes + iMisfitType).getValue(shotInd)) / misfitStorageL2.at(misfitStorageL2Size - 2 * numMisfitTypes + iMisfitType).getValue(shotInd);
            }
            auto misfitRatioMax = std::max_element(std::begin(misfitRatio), std::end(misfitRatio));
            misfitRatioMaxInd = std::distance(std::begin(misfitRatio), misfitRatioMax);
        }
        misfitTypeShots.setValue(shotInd, uniqueMisfitTypes.at(misfitRatioMaxInd));
        misfitSum = misfitStorageL2.at(misfitStorageL2Size - numMisfitTypes + misfitRatioMaxInd).getValue(shotInd)*count*misfitTypeShots.size();
    } else if ((misfitType.length() > 2 && misfitType.substr(numMisfitTypes+1, 1).compare("1") == 0) || saveMultiMisfits) {
        // for single misfitType in main code
        for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) {   
            if (uniqueMisfitTypes.at(iMisfitType) == misfitTypeShots.getValue(shotInd)) {
                misfitSum = misfitStorageL2.at(misfitStorageL2Size - numMisfitTypes + iMisfitType).getValue(shotInd)*count*misfitTypeShots.size();
                break;
            }
        }
    } else {
        // for single misfitType in steplength search
        for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) {   
            if (uniqueMisfitTypes.at(iMisfitType) == misfitTypeShots.getValue(shotInd)) {
                misfitSum *= misfitSum0Ratio.at(iMisfitType).getValue(shotInd);
                break;
            }
        }
    }
    
    return misfitSum/count/misfitTypeShots.size(); 
}

/*! \brief Calculate the adjoint sources
 *
 *
 \param adjointSources Receiver object which stores the adjoint sources (in- and output)
 \param receiversSyn Receiver object which stores the synthetic data 
 \param receiversObs Receiver object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSources(KITGPI::Acquisition::Receivers<ValueType> &adjointSources, KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs, scai::IndexType shotInd)
{      
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObs;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramAdj;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerSyn;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerObs;
    
    seismoHandlerSyn = receiversSyn.getSeismogramHandler();
    seismoHandlerObs = receiversObs.getSeismogramHandler();
    
    scai::IndexType misfitTypeShotL2 = misfitTypeShots.getValue(shotInd);
    for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
        seismogramObs = seismoHandlerObs.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
        switch (misfitTypeShotL2) {
        case 5:
            this->calcAdjointSeismogramL2EnvelopeWeighted(seismogramAdj, seismogramSyn, seismogramObs);
            break;
        case 6:
            this->calcAdjointSeismogramL2AGC(seismogramAdj, seismogramSyn, seismogramObs);
            break;
        case 7:
            this->calcAdjointSeismogramL2Normalized(seismogramAdj, seismogramSyn, seismogramObs);
            break;
        case 8:
            this->calcAdjointSeismogramL2Envelope(seismogramAdj, seismogramSyn, seismogramObs);
            break;
        case 9:
            this->calcAdjointSeismogramL2InstantaneousPhase(seismogramAdj, seismogramSyn, seismogramObs);
            break;
        default:
            this->calcAdjointSeismogramL2(seismogramAdj, seismogramSyn, seismogramObs);
            break;
        } 
        adjointSources.getSeismogramHandler().getSeismogram(seismogramAdj.getTraceType()) = seismogramAdj;
    } 
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;    
    if (seismogramSyn.getData().getNumRows()!=0) {
        // seismogramSyntemp *= seismogramObstemp.getData().maxNorm() / seismogramSyntemp.getData().maxNorm();
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows();  
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;    
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    if (seismogramSyn.getData().getNumRows()!=0) {
        seismogramAdj = seismogramSyn;
//         seismogramAdj *= seismogramObs.getData().maxNorm() / seismogramSyn.getData().maxNorm();
        seismogramAdj -= seismogramObs; 
        
        bool isSeismic = seismogramSyn.getIsSeismic();
        if ((isSeismic && seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) || (!isSeismic && seismogramSyn.getTraceTypeEM() != Acquisition::SeismogramTypeEM::HZ)) {
            // Seismic adjoint source P and EM adjoint source E times -1 for different reason. Seismic adjoint source P times -1 in an anti self-adjoint state equation, while EM adjoint source E times -1 in a self-adjoint state equation. Please see details in the manual.
            seismogramAdj *= -1;        
        }
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2EnvelopeWeighted(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramSyntemp.normalizeTrace(4);
        seismogramObstemp.normalizeTrace(4);
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows(); 
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;   
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2EnvelopeWeighted(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        // seismogramSyntemp *= seismogramObstemp.getData().maxNorm() / seismogramSyntemp.getData().maxNorm();
        seismogramAdj = seismogramSyntemp;
        scai::lama::DenseMatrix<ValueType> tempDataSyn;
        scai::lama::DenseMatrix<ValueType> tempDataSynEnvelope = seismogramSyntemp.getData();
        scai::lama::DenseMatrix<ValueType> tempDataSynRatio;
        scai::lama::DenseVector<ValueType> dataTrace; 
        ValueType waterLevel = 1e-3;
        Common::calcEnvelope(tempDataSynEnvelope);
        seismogramSyntemp.normalizeTrace(4);
        seismogramObstemp.normalizeTrace(4);
        tempDataSynRatio.binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::SUB, seismogramObstemp.getData());
        for (int i=0; i<tempDataSynEnvelope.getNumRows(); i++) {
            tempDataSynEnvelope.getRow(dataTrace, i);   
            if (dataTrace.maxNorm() != 0) {
                dataTrace += waterLevel * dataTrace.maxNorm();
            } else {
                dataTrace += waterLevel * waterLevel;
            }
            tempDataSynEnvelope.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
        }
        tempDataSynRatio.binaryOp(tempDataSynRatio, scai::common::BinaryOp::DIVIDE, tempDataSynEnvelope);
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, seismogramSyntemp.getData());
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, tempDataSynRatio);
        seismogramSyntemp.getData().binaryOp(tempDataSynRatio, scai::common::BinaryOp::SUB, seismogramSyntemp.getData());
        
        tempDataSynEnvelope.binaryOp(tempDataSynEnvelope, scai::common::BinaryOp::MULT, tempDataSynEnvelope);
        tempDataSyn = seismogramAdj.getData();
        Hilbert::HilbertFFT<ValueType> hilbertHandler;
        IndexType hLength = Common::calcNextPowTwo<ValueType>(tempDataSyn.getNumColumns());  
        hilbertHandler.setCoefficientLength(hLength);
        hilbertHandler.calcHilbertCoefficient();
        hilbertHandler.hilbert(tempDataSyn);
        tempDataSyn.binaryOp(tempDataSyn, scai::common::BinaryOp::DIVIDE, tempDataSynEnvelope);
        tempDataSyn.binaryOp(tempDataSyn, scai::common::BinaryOp::MULT, seismogramAdj.getData());
        tempDataSyn.binaryOp(tempDataSyn, scai::common::BinaryOp::MULT, tempDataSynRatio);
        hilbertHandler.hilbert(tempDataSyn);
        seismogramObstemp.getData() = tempDataSyn;   
        
        seismogramAdj = seismogramSyntemp - seismogramObstemp;  
        
        bool isSeismic = seismogramSyn.getIsSeismic();
        if ((isSeismic && seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) || (!isSeismic && seismogramSyn.getTraceTypeEM() != Acquisition::SeismogramTypeEM::HZ)) {
            seismogramAdj *= -1;        
        }
    }
}

/*! \brief Return the AGC weighted L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2AGC(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;    
    if (seismogramSyn.getData().getNumRows()!=0) {
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, seismogramSyn.getInverseAGC());   
        seismogramObstemp.getData().binaryOp(seismogramObstemp.getData(), scai::common::BinaryOp::MULT, seismogramObs.getInverseAGC());    
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows();  
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;    
}

/*! \brief Calculate the AGC weighted adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2AGC(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;    
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        // seismogramSyntemp *= seismogramObstemp.getData().maxNorm() / seismogramSyntemp.getData().maxNorm();
        seismogramAdj = seismogramSyntemp;
        seismogramAdj *= seismogramObstemp;
        scai::lama::DenseMatrix<ValueType> sumSynMultObs = seismogramAdj.getAGCSum();
        sumSynMultObs.binaryOp(sumSynMultObs, scai::common::BinaryOp::MULT, seismogramSyn.getInverseAGC()); 
        sumSynMultObs.binaryOp(sumSynMultObs, scai::common::BinaryOp::MULT, seismogramObs.getInverseAGC()); 
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, seismogramSyn.getInverseAGC());   
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, sumSynMultObs);       
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, seismogramSyn.getInverseAGC());  
        seismogramObstemp.getData().binaryOp(seismogramObstemp.getData(), scai::common::BinaryOp::MULT, seismogramObs.getInverseAGC());     
        seismogramObstemp.getData().binaryOp(seismogramObstemp.getData(), scai::common::BinaryOp::MULT, seismogramSyn.getInverseAGC());  
        
        seismogramAdj = seismogramSyntemp - seismogramObstemp; 
            
        bool isSeismic = seismogramSyn.getIsSeismic();
        if ((isSeismic && seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) || (!isSeismic && seismogramSyn.getTraceTypeEM() != Acquisition::SeismogramTypeEM::HZ)) {
            seismogramAdj *= -1;        
        }
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2Normalized(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramSyntemp.normalizeTrace(2);
        seismogramObstemp.normalizeTrace(2);
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows(); 
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;   
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2Normalized(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
  /*      // Groos 2013 and IFOS2D
        seismogramAdj = seismogramSyn;
        seismogramAdj *= seismogramObs;
        scai::lama::DenseVector<ValueType> tempL2NormSyn = seismogramSyntemp.getTraceL2norm();
        scai::lama::DenseVector<ValueType> tempL2NormObs = seismogramObstemp.getTraceL2norm();
        scai::lama::DenseVector<ValueType> sumSynMultObs = seismogramAdj.getTraceSum();
        tempL2NormSyn = 1 / tempL2NormSyn;
        tempL2NormObs = 1 / tempL2NormObs;     
        tempL2NormObs *= tempL2NormSyn;
        tempL2NormObs *= tempL2NormSyn;
        tempL2NormObs *= sumSynMultObs;
        seismogramSyntemp.normalizeTrace(2);
        seismogramObstemp.normalizeTrace(2);
        seismogramSyntemp.getData().scaleRows(tempL2NormObs);
        seismogramObstemp.getData().scaleRows(tempL2NormSyn);  */  
        // Choi 2012
        // seismogramSyntemp *= seismogramObstemp.getData().maxNorm() / seismogramSyntemp.getData().maxNorm();
        scai::lama::DenseVector<ValueType> tempL2NormSyn = seismogramSyntemp.getTraceL2norm();        
        Common::searchAndReplace<ValueType>(tempL2NormSyn, 0.0, 1.0, 5);
        tempL2NormSyn = 1 / tempL2NormSyn;
        seismogramSyntemp.normalizeTrace(2);
        seismogramObstemp.normalizeTrace(2);
        seismogramAdj = seismogramSyntemp;
        seismogramAdj *= seismogramObstemp;
        scai::lama::DenseVector<ValueType> sumSynMultObs = seismogramAdj.getTraceSum();
        seismogramSyntemp.getData().scaleRows(sumSynMultObs);       
        seismogramSyntemp.getData().scaleRows(tempL2NormSyn);       
        seismogramObstemp.getData().scaleRows(tempL2NormSyn);
        
        seismogramAdj = seismogramSyntemp - seismogramObstemp;       
        
        bool isSeismic = seismogramSyn.getIsSeismic();
        if ((isSeismic && seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) || (!isSeismic && seismogramSyn.getTraceTypeEM() != Acquisition::SeismogramTypeEM::HZ)) {
            seismogramAdj *= -1;        
        }
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2Envelope(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        // seismogramSyntemp *= seismogramObstemp.getData().maxNorm() / seismogramSyntemp.getData().maxNorm();
        Common::calcEnvelope(seismogramSyntemp.getData());
        Common::calcEnvelope(seismogramObstemp.getData());
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows(); 
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;   
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2Envelope(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        // seismogramSyntemp *= seismogramObstemp.getData().maxNorm() / seismogramSyntemp.getData().maxNorm();
        seismogramAdj = seismogramSyntemp;
        scai::lama::DenseMatrix<ValueType> tempDataSyn = seismogramSyntemp.getData();
        scai::lama::DenseMatrix<ValueType> tempDataObs = seismogramObstemp.getData();
        scai::lama::DenseVector<ValueType> dataTrace; 
        ValueType waterLevel = 1e-3;
        Common::calcEnvelope(tempDataSyn);
        Common::calcEnvelope(tempDataObs);
        tempDataObs.binaryOp(tempDataSyn, scai::common::BinaryOp::SUB, tempDataObs);
        for (int i=0; i<tempDataSyn.getNumRows(); i++) {
            tempDataSyn.getRow(dataTrace, i);   
            if (dataTrace.maxNorm() != 0) {
                dataTrace += waterLevel * dataTrace.maxNorm();
            } else {
                dataTrace += waterLevel * waterLevel;
            }
            tempDataSyn.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
        }
        tempDataObs.binaryOp(tempDataObs, scai::common::BinaryOp::DIVIDE, tempDataSyn);
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, tempDataObs);
        
        tempDataSyn = seismogramAdj.getData();    
        Hilbert::HilbertFFT<ValueType> hilbertHandler;
        IndexType hLength = Common::calcNextPowTwo<ValueType>(tempDataSyn.getNumColumns());  
        hilbertHandler.setCoefficientLength(hLength);
        hilbertHandler.calcHilbertCoefficient();
        hilbertHandler.hilbert(tempDataSyn);
        tempDataSyn.binaryOp(tempDataSyn, scai::common::BinaryOp::MULT, tempDataObs);
        hilbertHandler.hilbert(tempDataSyn);
        seismogramObstemp.getData() = tempDataSyn;  
        
        seismogramAdj = seismogramSyntemp - seismogramObstemp;      
        
        bool isSeismic = seismogramSyn.getIsSeismic();
        if ((isSeismic && seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) || (!isSeismic && seismogramSyn.getTraceTypeEM() != Acquisition::SeismogramTypeEM::HZ)) {
            seismogramAdj *= -1;        
        }
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2InstantaneousPhase(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        // seismogramSyntemp *= seismogramObstemp.getData().maxNorm() / seismogramSyntemp.getData().maxNorm();
        Common::calcInstantaneousPhaseResidual(seismogramSyntemp.getData(), seismogramObstemp.getData(), seismogramSyntemp.getData());
//         scai::IndexType phaseType = 1;
//         Common::calcInstantaneousPhase(seismogramSyntemp.getData(), phaseType);
//         Common::calcInstantaneousPhase(seismogramObstemp.getData(), phaseType);
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows(); 
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;   
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2InstantaneousPhase(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        // seismogramSyntemp *= seismogramObstemp.getData().maxNorm() / seismogramSyntemp.getData().maxNorm();
        seismogramAdj = seismogramSyntemp;
        scai::lama::DenseMatrix<ValueType> envelopeSyn = seismogramSyntemp.getData();
        scai::lama::DenseMatrix<ValueType> instantaneousPhaseSyn = seismogramSyntemp.getData();
        scai::lama::DenseMatrix<ValueType> instantaneousPhaseObs = seismogramObstemp.getData();
        scai::lama::DenseVector<ValueType> dataTrace; 
        ValueType waterLevel = 1e-3;
        Common::calcEnvelope(envelopeSyn);
        envelopeSyn.binaryOp(envelopeSyn, scai::common::BinaryOp::MULT, envelopeSyn);
        for (int i=0; i<envelopeSyn.getNumRows(); i++) {
            envelopeSyn.getRow(dataTrace, i);   
            if (dataTrace.maxNorm() != 0) {
                dataTrace += waterLevel * dataTrace.maxNorm();
            } else {
                dataTrace += waterLevel * waterLevel;
            }
            envelopeSyn.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
        }
//         scai::IndexType phaseType = 1;
//         Common::calcInstantaneousPhase(instantaneousPhaseObs, phaseType);
//         Common::calcInstantaneousPhase(instantaneousPhaseSyn, phaseType);
//         instantaneousPhaseObs.binaryOp(instantaneousPhaseObs, scai::common::BinaryOp::SUB, instantaneousPhaseSyn);
        Common::calcInstantaneousPhaseResidual(instantaneousPhaseObs, instantaneousPhaseObs, instantaneousPhaseSyn);
        
        instantaneousPhaseObs.binaryOp(instantaneousPhaseObs, scai::common::BinaryOp::DIVIDE, envelopeSyn);
        
        scai::lama::DenseMatrix<ValueType> tempDataSyn = seismogramSyntemp.getData();    
        Hilbert::HilbertFFT<ValueType> hilbertHandler;
        IndexType hLength = Common::calcNextPowTwo<ValueType>(tempDataSyn.getNumColumns());  
        hilbertHandler.setCoefficientLength(hLength);
        hilbertHandler.calcHilbertCoefficient();
        hilbertHandler.hilbert(tempDataSyn);
        seismogramSyntemp.getData().binaryOp(tempDataSyn, scai::common::BinaryOp::MULT, instantaneousPhaseObs);
        
        tempDataSyn = seismogramAdj.getData();
        tempDataSyn.binaryOp(tempDataSyn, scai::common::BinaryOp::MULT, instantaneousPhaseObs);
        hilbertHandler.hilbert(tempDataSyn);
        seismogramObstemp.getData() = tempDataSyn; 
        
        seismogramAdj = seismogramSyntemp + seismogramObstemp;
    
        bool isSeismic = seismogramSyn.getIsSeismic();
        if ((isSeismic && seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) || (!isSeismic && seismogramSyn.getTraceTypeEM() != Acquisition::SeismogramTypeEM::HZ)) {
            seismogramAdj *= -1;        
        }
    }
}

template class KITGPI::Misfit::MisfitL2<double>;
template class KITGPI::Misfit::MisfitL2<float>;
