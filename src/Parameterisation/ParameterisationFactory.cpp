#include "ParameterisationFactory.hpp"

template <typename ValueType>
typename KITGPI::Parameterisation::Parameterisation<ValueType>::ParameterisationPtr KITGPI::Parameterisation::Factory<ValueType>::Create(std::string type)
{

    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("visco") == 0, "Unkown type");

    if (type.compare("acoustic") == 0) {
        return ParameterisationPtr(new Acoustic<ValueType>);
    }
    if (type.compare("elastic") == 0) {
        return ParameterisationPtr(new Elastic<ValueType>);
    }
    if (type.compare("visco") == 0) {
        return ParameterisationPtr(new Viscoelastic<ValueType>);
    }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return ParameterisationPtr();
}

template class KITGPI::Parameterisation::Factory<float>;
template class KITGPI::Parameterisation::Factory<double>;
