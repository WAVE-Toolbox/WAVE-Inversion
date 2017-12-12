#include "Misfit.hpp"

template <typename ValueType>
scai::lama::Scalar Misfit<ValueType>::getMisfitSum(int iteration){
    return this->misfitShot.at(iteration).sum();
}

template <typename ValueType>
scai::lama::Scalar Misfit<ValueType>::getMisfitShot(int iteration, int shotNumber){
    return this->misfitShot.at(iteration).getValue(shotNumber);
}

template <typename ValueType>
void Misfit<ValueType>::add(scai::lama::DenseVector<ValueType> vector){
    this->misfitShot.push_back(vector);
}


void l1(){}
void l2(){}

template class Misfit<double>;
template class Misfit<float>;

