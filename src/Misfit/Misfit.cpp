#include "Misfit.hpp"

/*! \brief Set misfit
 *
 \param type misfitType 
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::init(std::string type, scai::IndexType numshots)
{    
    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    if (type.length() == 2) {
        for (int i = 0; i < numshots; i++)
        {
            misfitTypeHistory.push_back(type);
        }
    } else if (type.compare("l278") == 0) {
        scai::IndexType randMisfitTypeInd;
        std::vector<std::string> randMisfitTypes{"l2","l7","l8"};
        std::srand((int)time(0));
        for (int i = 0; i < numshots; i++) {
            randMisfitTypeInd = std::rand() % randMisfitTypes.size();
            misfitTypeHistory.push_back(randMisfitTypes[randMisfitTypeInd]);
        }
    }
}

/*! \brief get misfit
 *
 \param type misfitType 
 */
template <typename ValueType>
std::string KITGPI::Misfit::Misfit<ValueType>::getMisfitType(scai::IndexType shotInd)
{    
    return misfitTypeHistory.at(shotInd);
}

/*! \brief get misfit history
 *
 \param type misfitType 
 */
template <typename ValueType>
std::vector<std::string>  KITGPI::Misfit::Misfit<ValueType>::getMisfitTypeHistory()
{    
    return misfitTypeHistory;
}

/*! \brief set misfit history vector
 *
 \param type misfitType 
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::setMisfitTypeHistory(std::vector<std::string> setMisfitTypeHistory)
{    
    misfitTypeHistory = setMisfitTypeHistory;
}

/*! \brief Write to misfitType-file
*
\param comm Communicator
\param misfitTypeFilename Name of misfitType-file
\param uniqueShotNos unique Shot numbers
\param uniqueShotNosRand unique Shot numbers randomly
\param stage inversion stage
\param iteration inversion iteration
\param misfitType misfitType
*/
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::writeMisfitTypeToFile(scai::dmemo::CommunicatorPtr comm, std::string misfitTypeFilename, std::vector<scai::IndexType> uniqueShotNos, std::vector<scai::IndexType> uniqueShotNosRand, scai::IndexType stage, scai::IndexType iteration, std::string misfitType)
{      
    int myRank = comm->getRank();  
    if (misfitType.length() > 2 && myRank == MASTERGPI) {
        std::ofstream outputFile; 
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
        scai::IndexType shotIndTrue = 0;
        scai::IndexType shotNumber;
        for (unsigned i = 0; i < uniqueShotNosRand.size(); i++) { 
            shotNumber = uniqueShotNosRand[i];
            Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNos, shotNumber);
            if (i == 0) {
                outputFile << std::setw(9) << this->getMisfitType(shotIndTrue);
            } else {
                outputFile << std::setw(4) << this->getMisfitType(shotIndTrue);
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
        
    std::cout << "calcStablizingFunctional stablizingFunctionalType = " << stablizingFunctionalType << std::endl;
    std::cout << "calcStablizingFunctional stablizingFunctional = " << tempValue << std::endl;
    
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
