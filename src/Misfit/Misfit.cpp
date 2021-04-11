#include "Misfit.hpp"

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
 \param shotNumber Integer value which specifies the shot number
 */
template <typename ValueType>
ValueType KITGPI::Misfit::Misfit<ValueType>::getMisfitShot(int iteration, int shotNumber)
{
    return this->misfitStorage.at(iteration).getValue(shotNumber);
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

template class KITGPI::Misfit::Misfit<double>;
template class KITGPI::Misfit::Misfit<float>;
