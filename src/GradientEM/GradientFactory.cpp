#include "GradientFactory.hpp"

template <typename ValueType>
typename KITGPI::Gradient::GradientEM<ValueType>::GradientPtr KITGPI::Gradient::FactoryEM<ValueType>::Create(std::string type)
{
    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("viscoelastic") == 0 || type.compare("sh") == 0 || type.compare("emem") == 0 || type.compare("tmem") == 0 || type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0, "Unkown type");
    
    if (type.compare("emem") == 0 || type.compare("tmem") == 0) {
        return GradientPtr(new EMEM<ValueType>);
    }
    if (type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0) {
        return GradientPtr(new ViscoEMEM<ValueType>);
    }

    return GradientPtr();
}

template class KITGPI::Gradient::FactoryEM<float>;
template class KITGPI::Gradient::FactoryEM<double>;