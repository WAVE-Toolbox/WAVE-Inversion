#include "ZeroLagXcorrFactory.hpp"

using namespace scai;

template <typename ValueType>
typename KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::ZeroLagXcorrPtr KITGPI::ZeroLagXcorr::Factory<ValueType>::Create(std::string dimension, std::string type)
{

    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(dimension.compare("2d") == 0 || dimension.compare("3d") == 0, "Unkown dimension");
    SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("visco") == 0, "Unkown type");

    // 2D
    if (dimension.compare("2d") == 0 && type.compare("acoustic") == 0) {
        return ZeroLagXcorrPtr(new ZeroLagXcorr2Dacoustic<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("elastic") == 0) {
        return ZeroLagXcorrPtr(new ZeroLagXcorr2Delastic<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("visco") == 0) {
        //         return ZeroLagXcorrPtr(new ZeroLagXcorr2Dvisco<ValueType>);
        COMMON_THROWEXCEPTION("3Dviscoelastic convolution is not implemented yet.");
    }

    // 3D
    if (dimension.compare("3d") == 0 && type.compare("acoustic") == 0) {
        return ZeroLagXcorrPtr(new ZeroLagXcorr3Dacoustic<ValueType>);
    }
    if (dimension.compare("3d") == 0 && type.compare("elastic") == 0) {
        return ZeroLagXcorrPtr(new ZeroLagXcorr3Delastic<ValueType>);
    }
    if (dimension.compare("3d") == 0 && type.compare("visco") == 0) {
        //         return ZeroLagXcorrPtr(new ZeroLagXcorr3Dvisco<ValueType>);
        COMMON_THROWEXCEPTION("3Dviscoelastic convolution is not implemented yet.");
    }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return ZeroLagXcorrPtr();
}

template class KITGPI::ZeroLagXcorr::Factory<double>;
template class KITGPI::ZeroLagXcorr::Factory<float>;
