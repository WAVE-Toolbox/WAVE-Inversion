#include "Misfit.hpp"

/*! \brief Set misfit
 *
 \param type misfitType 
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::init(KITGPI::Configuration::Configuration config, std::vector<scai::IndexType> misfitTypeHistory, scai::IndexType numshots)
{    
    misfitTypeShots.clear();
    // transform to lower cases
    std::string misfitType = config.get<std::string>("misfitType");
    std::transform(misfitType.begin(), misfitType.end(), misfitType.begin(), ::tolower);
    if (misfitType.length() == 2 && config.getAndCatch("useRandomSource", 0) == 0) {
        for (int shotInd = 0; shotInd < numshots; shotInd++) {
            misfitTypeShots.push_back(misfitType);
        }
    } else if (misfitType.compare("l278") == 0) {
        scai::IndexType randMisfitTypeInd;
        scai::IndexType numshotsPerIt = numshots;
        scai::IndexType maxiterations = config.get<scai::IndexType>("maxIterations");
        std::vector<std::string> uniqueMisfitTypes{"l2","l7","l8"};
        std::srand((int)time(0));
        if (config.getAndCatch("useRandomSource", 0) != 0) {
            numshotsPerIt = config.get<scai::IndexType>("NumShotDomains"); 
        }
        scai::IndexType maxcount = maxiterations * numshotsPerIt / uniqueMisfitTypes.size() * 1.2;
        for (int shotInd = 0; shotInd < numshotsPerIt; shotInd++) {         
            randMisfitTypeInd = std::rand() % uniqueMisfitTypes.size();
            if (misfitTypeHistory[randMisfitTypeInd] >= maxcount) {
                shotInd--;
            } else {
                misfitTypeShots.push_back(uniqueMisfitTypes[randMisfitTypeInd]);
                misfitTypeHistory[randMisfitTypeInd]++;
            }
        }
    }
}

/*! \brief get misfit
 *
 \param type misfitType 
 */
template <typename ValueType>
std::string KITGPI::Misfit::Misfit<ValueType>::getMisfitTypeShot(scai::IndexType shotInd)
{    
    return misfitTypeShots.at(shotInd);
}

/*! \brief get misfit history
 *
 \param type misfitType 
 */
template <typename ValueType>
std::vector<std::string>  KITGPI::Misfit::Misfit<ValueType>::getMisfitTypeShots()
{    
    return misfitTypeShots;
}

/*! \brief set misfit history vector
 *
 \param type misfitType 
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::setMisfitTypeShots(std::vector<std::string> setMisfitTypeShots)
{    
    misfitTypeShots = setMisfitTypeShots;
}

/*! \brief Write to misfitType-file
*
\param comm Communicator
\param logFilename Name of log-file
\param stage inversion stage
\param iteration inversion iteration
\param misfitType misfitType
*/
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::writeMisfitTypeToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, scai::IndexType stage, scai::IndexType iteration, std::string misfitType)
{      
    int myRank = comm->getRank();  
    if (misfitType.length() > 2 && myRank == MASTERGPI) {
        std::ofstream outputFile; 
        std::string misfitTypeFilename = logFilename.substr(0, logFilename.length()-4);
        misfitTypeFilename += ".misfitType";
        misfitTypeFilename += logFilename.substr(logFilename.length()-4, 4);
        if (stage == 1 && iteration == 1) {
            outputFile.open(misfitTypeFilename);
            outputFile << "# MisfitType records during inversion\n"; 
            outputFile << "# random misfit type = " << misfitType << "\n"; 
            outputFile << "# Stage | Iteration | misfitTypes\n"; 
        } else {                    
            outputFile.open(misfitTypeFilename, std::ios_base::app);
            outputFile << std::scientific;
        }
        outputFile << std::setw(5) << stage << std::setw(10) << iteration;
        for (unsigned i = 0; i < misfitTypeShots.size(); i++) { 
            if (i == 0) {
                outputFile << std::setw(9) << misfitTypeShots[i];
            } else {
                outputFile << std::setw(4) << misfitTypeShots[i];
            }
        }
        outputFile << "\n";
        outputFile.close();
    }
}

/*! \brief Return the misfit summed over all shots. 
 * 
 * 
 \param iteration Integer value which specifies the iteration number
 */
template <typename ValueType>
ValueType KITGPI::Misfit::Misfit<ValueType>::getMisfitSum(int iteration)
{
    return this->misfitStorage.at(iteration).sum();
}

/*! \brief Return a vector which stores the misfit of each shot for one iteration
 *
 *
 \param iteration Integer value which specifies the iteration number
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Misfit::Misfit<ValueType>::getMisfitIt(int iteration)
{
    return this->misfitStorage.at(iteration);
}

/*! \brief Return the misfit of one shot
 *
 *
 \param iteration Integer value which specifies the iteration number
 \param shotInd Integer value which specifies the shot number
 */
template <typename ValueType>
ValueType KITGPI::Misfit::Misfit<ValueType>::getMisfitShot(int iteration, int shotInd)
{
    return this->misfitStorage.at(iteration).getValue(shotInd);
}

/*! \brief Add the misfit of one iteration to the misfit storage
 *
 *
 \param vector Vector which stores the misfit of all shots separately during one iteration
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::addToStorage(scai::lama::DenseVector<ValueType> vector)
{
    this->misfitStorage.push_back(vector);
}

/*! \brief Clear the misfit storage 
 *
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::clearStorage()
{
    this->misfitStorage.clear();
}

/*! \brief Set number of Relaxation Mechanisms
 */
template <typename ValueType>
ValueType KITGPI::Misfit::Misfit<ValueType>::calcStablizingFunctionalPerModel(scai::lama::DenseVector<ValueType> modelResidualVec, ValueType focusingParameter, int stablizingFunctionalType)
{
    scai::lama::DenseVector<ValueType> temp;
    scai::lama::DenseVector<ValueType> temp1;
    ValueType tempValue;
    tempValue = focusingParameter;
    tempValue *= modelResidualVec.maxNorm();
    tempValue = pow(tempValue, 2);  
            
    switch (stablizingFunctionalType) {
        case 1:
            // case 1: arctangent(AT) functional
            temp = scai::lama::pow(modelResidualVec, 2);
            temp += tempValue;
            temp = scai::lama::sqrt(temp);
            temp = scai::lama::atan(temp);
            
            break;
        case 2:
            // case 2: minimum support (MS) functional   
            temp = scai::lama::pow(modelResidualVec, 2);
            temp1 = temp + tempValue;
            temp /= temp1;
            
            break;
        default:
            // case 0: no stablizing functional
            temp = modelResidualVec;
            temp *= 0;
    }     
    tempValue = temp.l1Norm();
        
//     std::cout << "calcStablizingFunctional stablizingFunctionalType = " << stablizingFunctionalType << std::endl;
//     std::cout << "calcStablizingFunctional stablizingFunctional = " << tempValue << std::endl;
//     
    return tempValue; 
}

/*! \brief Get const reference to modelDerivativeX
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Misfit::Misfit<ValueType>::getModelDerivativeX()
{
    return (modelDerivativeX);
}

/*! \brief Set modelDerivativeX
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::setModelDerivativeX(scai::lama::Vector<ValueType> const &setModelDerivativeX)
{
    modelDerivativeX = setModelDerivativeX;
}

/*! \brief Get const reference to modelDerivativeY
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Misfit::Misfit<ValueType>::getModelDerivativeY()
{
    return (modelDerivativeY);
}

/*! \brief Set modelDerivativeY
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::setModelDerivativeY(scai::lama::Vector<ValueType> const &setModelDerivativeY)
{
    modelDerivativeY = setModelDerivativeY;
}

template class KITGPI::Misfit::Misfit<double>;
template class KITGPI::Misfit::Misfit<float>;
