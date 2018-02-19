#include "Misfit.hpp"

template <typename ValueType>
scai::lama::Scalar KITGPI::Misfit::Misfit<ValueType>::getMisfitSum(int iteration)
{
    return this->misfitStorage.at(iteration).sum();
}

template <typename ValueType>
scai::lama::Scalar KITGPI::Misfit::Misfit<ValueType>::getMisfitShot(int iteration, int shotNumber)
{
    return this->misfitStorage.at(iteration).getValue(shotNumber);
}

template <typename ValueType>
void KITGPI::Misfit::Misfit<ValueType>::addToStorage(scai::lama::DenseVector<ValueType> vector)
{
    this->misfitStorage.push_back(vector);
}

template class KITGPI::Misfit::Misfit<double>;
template class KITGPI::Misfit::Misfit<float>;
