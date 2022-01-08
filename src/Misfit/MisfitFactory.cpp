#include "MisfitFactory.hpp"

template <typename ValueType>
typename KITGPI::Misfit::Misfit<ValueType>::MisfitPtr KITGPI::Misfit::Factory<ValueType>::Create(std::string type)
{
    
    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(type.compare("l2") == 0 || type.compare("l3") == 0 || type.compare("l4") == 0 || type.compare("l5") == 0 || type.compare("l6") == 0 || type.compare("l7") == 0 || type.compare("l8") == 0 || type.compare("l9") == 0 || type.length() > 2, "Unknown type");

    if (type.compare("l2") == 0 || type.compare("l3") == 0 || type.compare("l4") == 0 || type.compare("l5") == 0 || type.compare("l6") == 0 || type.compare("l7") == 0 || type.compare("l8") == 0 || type.compare("l9") == 0 || type.length() > 2) {
        return MisfitPtr(new MisfitL2<ValueType>);
    }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return MisfitPtr();
}

template class KITGPI::Misfit::Factory<float>;
template class KITGPI::Misfit::Factory<double>;
