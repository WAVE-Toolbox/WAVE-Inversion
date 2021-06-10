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
    std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::tolower);
    scai::IndexType misfitTypeLength = misfitType.length() - 2;
    std::string temp1 = misfitType.substr(1, misfitType.length());
    scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();                 // default context, set by environment variable SCAI_CONTEXT 
    scai::lama::DenseVector<ValueType> temp(numshots, 0, ctx);  
    misfitTypeShots = temp;
    if (misfitTypeLength > 1) {
        scai::lama::DenseVector<ValueType> misfitPerIt(numshots, 0, ctx);
        for (int iMisfitType = 0; iMisfitType < misfitTypeLength; iMisfitType++) {  
            misfitStorageL2.push_back(misfitPerIt); 
        }
        std::vector<scai::IndexType> temp{2, 7, 8};
        uniqueMisfitTypes = temp;
    } else {
        scai::IndexType misfitTypeNo = atoi(temp1.c_str());
        misfitTypeShots = misfitTypeNo;
    }
    if (misfitType.compare("l2781") == 0) {
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
    if (misfitType.length() > 2 && myRank == MASTERGPI) {
        std::ofstream outputFile; 
        std::string misfitTypeFilename = logFilename.substr(0, logFilename.length()-4) + ".misfitType" + logFilename.substr(logFilename.length()-4, 4);
        if (stage == 1 && iteration == 1) {
            outputFile.open(misfitTypeFilename);
            outputFile << "# MisfitType records during inversion\n"; 
            outputFile << "# misfit type = " << misfitType << "\n"; 
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
    scai::IndexType misfitTypeLength = misfitType.length() - 2;
    int myRank = comm->getRank();  
    if (misfitTypeLength > 1 && myRank == MASTERGPI) {
        std::ofstream outputFile; 
        std::string misfitTypeFilename = logFilename.substr(0, logFilename.length()-3) + misfitType + logFilename.substr(logFilename.length()-4, 4);
        if (stage == 1 && iteration == 1) {
            outputFile.open(misfitTypeFilename);
            outputFile << "# Misfit records during inversion\n"; 
            outputFile << "# misfit type = " << misfitType << "\n"; 
            outputFile << "# Stage | Iteration"; 
            for (int iMisfitType = 0; iMisfitType < misfitTypeLength; iMisfitType++) { 
                outputFile << " | " << std::setw(7) << uniqueMisfitTypes.at(iMisfitType) << std::setw(7);
            }
            outputFile << "\n"; 
        } else {                    
            outputFile.open(misfitTypeFilename, std::ios_base::app);
            outputFile << std::scientific;
        }
        outputFile << std::setw(5) << stage << std::setw(10) << iteration;
        
        scai::IndexType misfitStorageL2Size = misfitStorageL2.size();
        for (int iMisfitType = 0; iMisfitType < misfitTypeLength; iMisfitType++) { 
            scai::lama::DenseVector<ValueType> misfitPerIt = misfitStorageL2.at(misfitStorageL2Size - misfitTypeLength + iMisfitType);
            if (iMisfitType == 0) {
                outputFile << std::setw(18) << misfitPerIt.sum();
            } else {
                outputFile << std::setw(14) << misfitPerIt.sum();
            }
        }
        outputFile << "\n";
        outputFile.close(); 
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
    scai::IndexType misfitTypeLength = misfitType.length() - 2;
    if (misfitTypeLength > 1) {
        for (int iMisfitType = 0; iMisfitType < misfitTypeLength; iMisfitType++) { 
            scai::IndexType misfitStorageL2Size = misfitStorageL2.size();
            if (commInterShot->getRank() == MASTERGPI) {
                std::cout << "\n misfitPerIt sumShotDomain = ";
                for (int shotInd = 0; shotInd < misfitTypeShots.size(); shotInd++) { 
                    std::cout << misfitStorageL2.at(misfitStorageL2Size - misfitTypeLength + iMisfitType).getValue(shotInd) << "\t";
                } 
                std::cout << std::endl;
            }
            commInterShot->sumArray(misfitStorageL2.at(misfitStorageL2Size - misfitTypeLength + iMisfitType).getLocalValues());   
            if (commInterShot->getRank() == MASTERGPI) {
                std::cout << "\n misfitPerIt sumShotDomain sumArray = ";
                for (int shotInd = 0; shotInd < misfitTypeShots.size(); shotInd++) { 
                    std::cout << misfitStorageL2.at(misfitStorageL2Size - misfitTypeLength + iMisfitType).getValue(shotInd) << "\t";
                }      
                std::cout << std::endl;
            }
        }
        if (misfitType.compare("l2782") == 0) {
            if (commInterShot->getRank() == MASTERGPI) {
                std::cout << "\n misfitTypeShots sumShotDomain = ";
                for (int shotInd = 0; shotInd < misfitTypeShots.size(); shotInd++) { 
                    std::cout << misfitTypeShots.getValue(shotInd) << "\t";
                } 
                std::cout << std::endl;
            }
            commInterShot->sumArray(misfitTypeShots.getLocalValues()); 
//             misfitTypeShots /= misfitTypeShots.size();
            if (commInterShot->getRank() == MASTERGPI) {
                std::cout << "\n misfitTypeShots sumShotDomain sumArry = ";
                for (int shotInd = 0; shotInd < misfitTypeShots.size(); shotInd++) { 
                    std::cout << misfitTypeShots.getValue(shotInd) << "\t";
                } 
                std::cout << std::endl;
            }
        }
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
    return 0; 
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
}

/*! \brief Return the L2-norm of seismograms stored in the given receiver objects (note: the misfit is summed over all components, i.e. vx,vy,vz,p)
 *
 *
 \param receiversSyn Receiver object which stores the synthetic data 
 \param receiversObs Receiver object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversSyn, KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversObs, scai::IndexType shotInd)
{            
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObs;
    KITGPI::Acquisition::SeismogramHandlerEM<ValueType> seismoHandlerSyn;
    KITGPI::Acquisition::SeismogramHandlerEM<ValueType> seismoHandlerObs;
    
    ValueType misfit = 0;
    ValueType misfitSum = 0;
    scai::IndexType count = 0;

    seismoHandlerSyn = receiversSyn.getSeismogramHandler();
    seismoHandlerObs = receiversObs.getSeismogramHandler();
          
    /* Note that the misfit of different components (p,vx,vy,vz) is summed up. If p and v? is used at the same time, this could cause problems because they could have different scales.
       For different velocity components it should be ok. */
    
    scai::IndexType misfitStorageL2Size = 0;
    scai::IndexType misfitTypeShotL2 = misfitTypeShots.getValue(shotInd);
    scai::IndexType misfitTypeLength = misfitType.length() - 2;
    if (misfitTypeLength == 0)
        misfitTypeLength += 1;
    for (int iMisfitType = 0; iMisfitType < misfitTypeLength; iMisfitType++) {  
        if (misfitTypeLength > 1) 
           misfitTypeShotL2 = uniqueMisfitTypes.at(iMisfitType);
            
        for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
            seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<Acquisition::SeismogramTypeEM>(i));
            seismogramObs = seismoHandlerObs.getSeismogram(static_cast<Acquisition::SeismogramTypeEM>(i));
            if (misfitTypeShotL2 == 2) {
                misfit = this->calcL2(seismogramSyn, seismogramObs);
            } else if (misfitTypeShotL2 == 7) {
                misfit = this->calcL2Normalized(seismogramSyn, seismogramObs);
            } else if (misfitTypeShotL2 == 8) {
                misfit = this->calcL2Envelope(seismogramSyn, seismogramObs);
            } 
            misfitSum += misfit;
            if (misfit != 0) count++;
        }  
        
        if (misfitTypeLength > 1) {
            misfitStorageL2Size = misfitStorageL2.size();
            misfitStorageL2.at(misfitStorageL2Size - misfitTypeLength + iMisfitType).setValue(shotInd, misfitSum/count/misfitTypeShots.size());
        }
        if (iMisfitType < misfitTypeLength - 1) {
            misfitSum = 0;
            count = 0;
        }
    }
    if (misfitType.compare("l2782") == 0) {
        scai::IndexType misfitRatioMaxInd = 0;
        if (misfitStorageL2Size > misfitTypeLength) { // iteration > 0
            std::vector<ValueType> misfitRatio(misfitTypeLength, 0);  
            for (int iMisfitType = 0; iMisfitType < misfitTypeLength; iMisfitType++) {           
                misfitRatio.at(iMisfitType) = (misfitStorageL2.at(misfitStorageL2Size - 2 * misfitTypeLength + iMisfitType).getValue(shotInd) - misfitStorageL2.at(misfitStorageL2Size - misfitTypeLength + iMisfitType).getValue(shotInd)) / misfitStorageL2.at(misfitStorageL2Size - 2 * misfitTypeLength + iMisfitType).getValue(shotInd);
            }
            auto misfitRatioMax = std::max_element(std::begin(misfitRatio), std::end(misfitRatio));
            misfitRatioMaxInd = std::distance(std::begin(misfitRatio), misfitRatioMax);
        }
        misfitTypeShots.setValue(shotInd, uniqueMisfitTypes.at(misfitRatioMaxInd));
        misfitSum = misfitStorageL2.at(misfitStorageL2Size - misfitTypeLength + misfitRatioMaxInd).getValue(shotInd)*count*misfitTypeShots.size();
    } else if (misfitType.compare("l2781") == 0) {
        for (int iMisfitType = 0; iMisfitType < misfitTypeLength; iMisfitType++) {   
            if (uniqueMisfitTypes.at(iMisfitType) == misfitTypeShots.getValue(shotInd)) {
                misfitSum = misfitStorageL2.at(misfitStorageL2Size - misfitTypeLength + iMisfitType).getValue(shotInd)*count*misfitTypeShots.size();
                break;
            }
        }
    }
    
    return misfitSum/count/misfitTypeShots.size(); 
}

/*! \brief Calculate the adjoint sources
 *
 *
 \param adjointSourcesEM Receiver object which stores the adjoint sources (in- and output)
 \param receiversSyn Receiver object which stores the synthetic data 
 \param receiversObs Receiver object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSources(KITGPI::Acquisition::ReceiversEM<ValueType> &adjointSourcesEM, KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversSyn, KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversObs, scai::IndexType shotInd)
{
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObs;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramAdj;
    KITGPI::Acquisition::SeismogramHandlerEM<ValueType> seismoHandlerSyn;
    KITGPI::Acquisition::SeismogramHandlerEM<ValueType> seismoHandlerObs;
    
    seismoHandlerSyn = receiversSyn.getSeismogramHandler();
    seismoHandlerObs = receiversObs.getSeismogramHandler();
    
    scai::IndexType misfitTypeShotL2 = misfitTypeShots.getValue(shotInd);
    for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<Acquisition::SeismogramTypeEM>(i));
        seismogramObs = seismoHandlerObs.getSeismogram(static_cast<Acquisition::SeismogramTypeEM>(i));
        if (misfitTypeShotL2 == 2) {
            this->calcAdjointSeismogramL2(seismogramAdj, seismogramSyn, seismogramObs);
        } else if (misfitTypeShotL2 == 7) {
            this->calcAdjointSeismogramL2Normalized(seismogramAdj, seismogramSyn, seismogramObs);
        } else if (misfitTypeShotL2 == 8) {
            this->calcAdjointSeismogramL2Envelope(seismogramAdj, seismogramSyn, seismogramObs);
        } 
        adjointSourcesEM.getSeismogramHandler().getSeismogram(seismogramAdj.getTraceType()) = seismogramAdj;
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

    seismogramAdj = seismogramSyn - seismogramObs; 
    
    if (seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) {
        seismogramAdj *= -1;        
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2(KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    
    if (seismogramSyn.getData().getNumRows()!=0) {
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
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2(KITGPI::Acquisition::SeismogramEM<ValueType> &seismogramAdj, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    seismogramAdj = seismogramSyn - seismogramObs; 
    
    if (seismogramSyn.getTraceType() != Acquisition::SeismogramTypeEM::HZ) {
        seismogramAdj *= -1;        
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

    // see Groos L. 2D full waveform inversion of shallow seismic Rayleigh waves[D]. Verlag nicht ermittelbar, 2013 for details
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramSyntemp.normalizeTraceL2();
        seismogramObstemp.normalizeTraceL2();
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

    // see Groos L. 2D full waveform inversion of shallow seismic Rayleigh waves[D]. Verlag nicht ermittelbar, 2013 for details
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramAdj = seismogramSyn;
        seismogramAdj *= seismogramObs;
        scai::lama::DenseVector<ValueType> tempL2NormSyn = seismogramSyntemp.getTraceL2norm();
        scai::lama::DenseVector<ValueType> tempL2NormObs = seismogramObstemp.getTraceL2norm();
        scai::lama::DenseVector<ValueType> tempSumSynObs = seismogramAdj.getTraceSum();
        tempL2NormSyn = 1 / tempL2NormSyn;
        tempL2NormObs = 1 / tempL2NormObs;     
        tempL2NormObs *= tempL2NormSyn;
        tempL2NormObs *= tempL2NormSyn;
        tempL2NormObs *= tempSumSynObs;
        seismogramSyntemp.normalizeTraceL2();
        seismogramObstemp.normalizeTraceL2();
        seismogramSyntemp.getData().scaleRows(tempL2NormObs);
        seismogramObstemp.getData().scaleRows(tempL2NormSyn);      
    }
    seismogramAdj = seismogramSyntemp - seismogramObstemp;
    
    if (seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) {
        seismogramAdj *= -1;        
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2Normalized(KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramSyntemp.normalizeTraceL2();
        seismogramObstemp.normalizeTraceL2();
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
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2Normalized(KITGPI::Acquisition::SeismogramEM<ValueType> &seismogramAdj, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramAdj = seismogramSyn;
        seismogramAdj *= seismogramObs;
        scai::lama::DenseVector<ValueType> tempL2NormSyn = seismogramSyntemp.getTraceL2norm();
        scai::lama::DenseVector<ValueType> tempL2NormObs = seismogramObstemp.getTraceL2norm();
        scai::lama::DenseVector<ValueType> tempSumSynObs = seismogramAdj.getTraceSum();
        tempL2NormSyn = 1 / tempL2NormSyn;
        tempL2NormObs = 1 / tempL2NormObs;     
        tempL2NormObs *= tempL2NormSyn;
        tempL2NormObs *= tempL2NormSyn;
        tempL2NormObs *= tempSumSynObs;
        seismogramSyntemp.normalizeTraceL2();
        seismogramObstemp.normalizeTraceL2();
        seismogramSyntemp.getData().scaleRows(tempL2NormObs);
        seismogramObstemp.getData().scaleRows(tempL2NormSyn);      
    }
    seismogramAdj = seismogramSyntemp - seismogramObstemp;       
    
    if (seismogramSyn.getTraceType() != Acquisition::SeismogramTypeEM::HZ) {
        seismogramAdj *= -1;        
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2Envelope(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    // see Yuan Y O, Simons F J, Bozdağ E. Multiscale adjoint waveform tomography for surface and body waves[J]. Geophysics, 2015, 80(5): R281-R302 for details
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        Common::envelope(seismogramSyntemp.getData());
        Common::envelope(seismogramObstemp.getData());
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

    // see Yuan Y O, Simons F J, Bozdağ E. Multiscale adjoint waveform tomography for surface and body waves[J]. Geophysics, 2015, 80(5): R281-R302 for details
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramAdj = seismogramSyn;
        scai::lama::DenseMatrix<ValueType> tempDataSyn = seismogramSyntemp.getData();
        scai::lama::DenseMatrix<ValueType> tempDataObs = seismogramObstemp.getData();
        scai::lama::DenseVector<ValueType> dataTrace; 
        ValueType waterLevel = 1.0;
        Common::envelope(tempDataSyn);
        Common::envelope(tempDataObs);
        for (int i=0; i<tempDataObs.getNumRows(); i++) {
            tempDataObs.getRow(dataTrace, i);  
            if (i==0) {
                waterLevel = dataTrace.maxNorm();                
            } else {
                if (waterLevel > dataTrace.maxNorm()) waterLevel = dataTrace.maxNorm(); 
            }
        }
        waterLevel *= 1e-1;
        tempDataObs.binaryOp(tempDataSyn, scai::common::BinaryOp::SUB, tempDataObs);
        for (int i=0; i<tempDataSyn.getNumRows(); i++) {
            tempDataSyn.getRow(dataTrace, i);   
            dataTrace += waterLevel;
            tempDataSyn.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
        }
        tempDataObs.binaryOp(tempDataObs, scai::common::BinaryOp::DIVIDE, tempDataSyn);
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, tempDataObs);
        
        tempDataSyn = seismogramAdj.getData();
        Common::hilbert(tempDataSyn);
        tempDataSyn.binaryOp(tempDataSyn, scai::common::BinaryOp::MULT, tempDataObs);
        Common::hilbert(tempDataSyn);
        seismogramObstemp.getData() = tempDataSyn;  
    }
    seismogramAdj = seismogramSyntemp - seismogramObstemp;
    
    if (seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) {
        seismogramAdj *= -1;        
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2Envelope(KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        Common::envelope(seismogramSyntemp.getData());
        Common::envelope(seismogramObstemp.getData());
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
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2Envelope(KITGPI::Acquisition::SeismogramEM<ValueType> &seismogramAdj, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramAdj = seismogramSyn;
        scai::lama::DenseMatrix<ValueType> tempDataSyn = seismogramSyntemp.getData();
        scai::lama::DenseMatrix<ValueType> tempDataObs = seismogramObstemp.getData();
        scai::lama::DenseVector<ValueType> dataTrace; 
        ValueType waterLevel = 1.0;
        Common::envelope(tempDataSyn);
        Common::envelope(tempDataObs);
        for (int i=0; i<tempDataObs.getNumRows(); i++) {
            tempDataObs.getRow(dataTrace, i);  
            if (i==0) {
                waterLevel = dataTrace.maxNorm();                
            } else {
                if (waterLevel > dataTrace.maxNorm()) waterLevel = dataTrace.maxNorm(); 
            }
        }
        waterLevel *= 1e-1;
        tempDataObs.binaryOp(tempDataSyn, scai::common::BinaryOp::SUB, tempDataObs);
        for (int i=0; i<tempDataSyn.getNumRows(); i++) {
            tempDataSyn.getRow(dataTrace, i);   
            dataTrace += waterLevel;
            tempDataSyn.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
        }
        tempDataObs.binaryOp(tempDataObs, scai::common::BinaryOp::DIVIDE, tempDataSyn);
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, tempDataObs);
        
        tempDataSyn = seismogramAdj.getData();
        Common::hilbert(tempDataSyn);
        tempDataSyn.binaryOp(tempDataSyn, scai::common::BinaryOp::MULT, tempDataObs);
        Common::hilbert(tempDataSyn);
        seismogramObstemp.getData() = tempDataSyn;  
    }
    seismogramAdj = seismogramSyntemp - seismogramObstemp;      
    
    if (seismogramSyn.getTraceType() != Acquisition::SeismogramTypeEM::HZ) {
        seismogramAdj *= -1;        
    }
}

template class KITGPI::Misfit::MisfitL2<double>;
template class KITGPI::Misfit::MisfitL2<float>;
