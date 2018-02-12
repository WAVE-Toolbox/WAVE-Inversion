#include "Misfit.hpp"

/*! \brief Get the misfit sum of all shots for a given iteration
 *
 *
 \param iteration Number of iteration
 */
template <typename ValueType>
scai::lama::Scalar Misfit<ValueType>::getMisfitSum(int iteration)
{
    return this->misfitShot.at(iteration).sum();
}

/*! \brief Get the misfit of a single shot in one iteration
 *
 *
 \param iteration Number of iteration
 \param shotNumber Number of shot
 */
template <typename ValueType>
scai::lama::Scalar Misfit<ValueType>::getMisfitShot(int iteration, int shotNumber)
{
    return this->misfitShot.at(iteration).getValue(shotNumber);
}

/*! \brief Add misfits
 *
 *
 \param vector Vector of misfits
 */
template <typename ValueType>
void Misfit<ValueType>::add(scai::lama::DenseVector<ValueType> vector)
{
    this->misfitShot.push_back(vector);
}

void l1() {}
void l2() {}

template class Misfit<double>;
template class Misfit<float>;
