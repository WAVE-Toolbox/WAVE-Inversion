#include "TaperFactory.hpp"

template <typename ValueType>
typename KITGPI::Taper::Taper<ValueType>::TaperPtr KITGPI::Taper::Factory<ValueType>::Create(std::string dimension)
{

    // transform to lower cases
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(dimension.compare("1d") == 0 || dimension.compare("2d") == 0, "Unkown dimension");

    // 1D
    if (dimension.compare("1d") == 0) {
        return TaperPtr(new Taper1D<ValueType>);
    }

    // 2D
    if (dimension.compare("2d") == 0) {
        return TaperPtr(new Taper2D<ValueType>);
    }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return TaperPtr();
};

template class KITGPI::Taper::Factory<double>;
template class KITGPI::Taper::Factory<float>;