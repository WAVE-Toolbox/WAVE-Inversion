#include "GradientFactory.hpp"

template <typename ValueType>
typename KITGPI::Gradient::Gradient<ValueType>::GradientPtr KITGPI::Gradient::Factory<ValueType>::Create(std::string type)
{
    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(type.compare("sh") == 0 || type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("visco") == 0, "Unkown type");

    if (type.compare("sh") == 0) {
        return GradientPtr(new SH<ValueType>);
    }
    if (type.compare("acoustic") == 0) {
        return GradientPtr(new Acoustic<ValueType>);
    }
    if (type.compare("elastic") == 0) {
        return GradientPtr(new Elastic<ValueType>);
    }
    if (type.compare("visco") == 0) {
        return GradientPtr(new Viscoelastic<ValueType>);
    }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return GradientPtr();
}

template class KITGPI::Gradient::Factory<float>;
template class KITGPI::Gradient::Factory<double>;
