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

template class KITGPI::Misfit::Misfit<double>;
template class KITGPI::Misfit::Misfit<float>;
