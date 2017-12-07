#include "Misfit.hpp"

template <typename ValueType>
scai::lama::Scalar Misfit<ValueType>::get(int element){
    return this->misfit.at(element);
}

template <typename ValueType>
void Misfit<ValueType>::set(int element, scai::lama::Scalar value){
    this->misfit.at(element) = value;
}

template <typename ValueType>
void Misfit<ValueType>::add(scai::lama::Scalar value){
    this->misfit.push_back(value);
}

void l1(){}
void l2(){}


template class Misfit<double>;
template class Misfit<float>;

