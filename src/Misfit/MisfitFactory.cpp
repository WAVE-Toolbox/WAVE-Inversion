#include "MisfitFactory.hpp"

template <typename ValueType>
typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr KITGPI::Misfit::Factory<ValueType>::Create(std::string type)
{
    
    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(type.compare("l2") == 0 || type.compare("l7") == 0 || type.compare("l8") == 0 || type.compare("l2781") == 0 || type.compare("l2782") == 0, "Unknown type");

    if (type.compare("l2") == 0 || type.compare("l7") == 0 || type.compare("l8") == 0 || type.compare("l2781") == 0 || type.compare("l2782") == 0) {
        return MisfitPtr(new MisfitL2<ValueType>);
    }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return MisfitPtr();
}

template class KITGPI::Misfit::Factory<float>;
template class KITGPI::Misfit::Factory<double>;
