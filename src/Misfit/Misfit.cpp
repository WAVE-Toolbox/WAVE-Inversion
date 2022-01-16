#include "Misfit.hpp"

/*! \brief get misfitType of all shots
 *
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Misfit::Misfit<ValueType>::getMisfitTypeShots()
{    
    return misfitTypeShots;
}

/*! \brief set misfitType of all shots
 *
 \param setMisfitTypeShots misfitType 
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::setMisfitTypeShots(scai::lama::DenseVector<ValueType> setMisfitTypeShots)
{    
    misfitTypeShots = setMisfitTypeShots;
}

/*! \brief get misfitSum0Ratio of all shots
 *
 */
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> KITGPI::Misfit::Misfit<ValueType>::getMisfitSum0Ratio()
{    
    return misfitSum0Ratio;
}

/*! \brief set misfitSum0Ratio of all shots
 *
 \param setMisfitTypeShots misfitType 
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::setMisfitSum0Ratio(std::vector<scai::lama::DenseVector<ValueType>> setMisfitSum0Ratio)
{    
    misfitSum0Ratio = setMisfitSum0Ratio;
}

/*! \brief Return the misfit residual vector over all common shots at iteration1 and iteration2. 
 \param iteration1 previous iteration
 \param iteration2 current iteration
 */
template <typename ValueType>
ValueType KITGPI::Misfit::Misfit<ValueType>::getMisfitResidualMax(int iteration1, int iteration2)
{
    ValueType misfitResidualMax = 0;    
    SCAI_ASSERT_ERROR(iteration1 < iteration2, "iteration1 >= iteration2");
    for (int iMisfitType = 0; iMisfitType < numMisfitTypes; iMisfitType++) {  
        scai::lama::DenseVector<ValueType> misfitPerIt1 = misfitStorageL2.at(iteration1 * numMisfitTypes + iMisfitType);
        scai::lama::DenseVector<ValueType> misfitPerIt2 = misfitStorageL2.at(iteration2 * numMisfitTypes + iMisfitType);
        scai::IndexType numshots = misfitPerIt1.size();
        ValueType misfit1 = 0;
        ValueType misfit2 = 0;
        for (int shotInd = 0; shotInd < numshots; shotInd++) {  
            if (misfitPerIt1.getValue(shotInd) != 0 && misfitPerIt2.getValue(shotInd) != 0) {
                misfit1 += misfitPerIt1.getValue(shotInd);
                misfit2 += misfitPerIt2.getValue(shotInd);
            }
        }
        SCAI_ASSERT_ERROR(misfit1 != 0, "misfit1 == 0 in the " + std::to_string(iMisfitType) + " misfitType");
        if (misfitResidualMax < (misfit1 - misfit2) / misfit1)
            misfitResidualMax = (misfit1 - misfit2) / misfit1;
    }
    return misfitResidualMax;
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

/*! \brief Add the crossGradientMisfit of one iteration to the crossGradientMisfit storage
 *
 \param crossGradientMisfit the crossGradientMisfit during one iteration
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::addToCrossGradientMisfitStorage(ValueType crossGradientMisfit)
{
    this->crossGradientMisfitStorage.push_back(crossGradientMisfit);
}

/*! \brief get the crossGradientMisfit at iteration
 *
 */
template <typename ValueType>
ValueType KITGPI::Misfit::Misfit<ValueType>::getCrossGradientMisfit(int iteration)
{
    return this->crossGradientMisfitStorage.at(iteration);
}

/*! \brief Clear the misfit storage 
 *
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::clearStorage()
{
    this->misfitStorage.clear();
    this->misfitStorageL2.clear();
    this->crossGradientMisfitStorage.clear();
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

/*! \brief Overloading = Operation
 *
 \param rhs Misfit which is copied.
 */
template <typename ValueType>
KITGPI::Misfit::Misfit<ValueType> &KITGPI::Misfit::Misfit<ValueType>::operator=(KITGPI::Misfit::Misfit<ValueType> const &rhs)
{
    misfitTypeShots = rhs.misfitTypeShots;
    misfitSum0Ratio = rhs.misfitSum0Ratio;
    uniqueMisfitTypes = rhs.uniqueMisfitTypes;
    modelDerivativeX = rhs.modelDerivativeX;
    modelDerivativeY = rhs.modelDerivativeY;
    fkHandler = rhs.fkHandler;
    nFFT = rhs.nFFT;
    
    return *this;
}

/*! \brief Calculate reflection sources
 *
 \param sourcesReflect sources for reflection.
 \param reflectivity reflectivity model.
 */
template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::calcReflectSources(KITGPI::Acquisition::Receivers<ValueType> &sourcesReflect, scai::lama::DenseVector<ValueType> reflectivity)
{
    for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        if (sourcesReflect.getSeismogramHandler().getSeismogram(static_cast<Acquisition::SeismogramType>(i)).getData().getNumRows()!=0) {
            if (static_cast<Acquisition::SeismogramTypeEM>(i) == Acquisition::SeismogramTypeEM::HZ) {
                reflectivity *= -1;        
            }
            sourcesReflect.getSeismogramHandler().getSeismogram(static_cast<Acquisition::SeismogramType>(i)) *= reflectivity;
        }
    }     
}

template class KITGPI::Misfit::Misfit<double>;
template class KITGPI::Misfit::Misfit<float>;
