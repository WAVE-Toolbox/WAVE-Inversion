#include "MisfitL2.hpp"

template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calc()
{
    return;   
}

// template <typename ValueType>
// void KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::Seismogram<ValueType> const &seismogram1, KITGPI::Acquisition::Seismogram<ValueType> const &seismogram2)
// {
//     return;   
// }

// template <typename ValueType>
// void KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::Seismogram<ValueType> const &seismogram1, KITGPI::Acquisition::Seismogram<ValueType> const &seismogram2)
// {
//     return;
//     
// }


template class KITGPI::Misfit::MisfitL2<double>;
template class KITGPI::Misfit::MisfitL2<float>;
