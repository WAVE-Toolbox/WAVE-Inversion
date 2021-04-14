#include "ZeroLagXcorrFactory.hpp"

using namespace scai;

template <typename ValueType>
typename KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::ZeroLagXcorrPtr KITGPI::ZeroLagXcorr::FactoryEM<ValueType>::Create(std::string dimensionEM, std::string type)
{

    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    std::transform(dimensionEM.begin(), dimensionEM.end(), dimensionEM.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(dimensionEM.compare("2d") == 0 || dimensionEM.compare("3d") == 0, "Unkown dimensionEM");
    SCAI_ASSERT_ERROR(type.compare("emem") == 0 || type.compare("tmem") == 0 || type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0, "Unkown type");

    // 2D
    if (dimensionEM.compare("2d") == 0 && type.compare("emem") == 0) {
        return ZeroLagXcorrPtr(new ZeroLagXcorr2Demem<ValueType>);
    }
    if (dimensionEM.compare("2d") == 0 && type.compare("tmem") == 0) {
        return ZeroLagXcorrPtr(new ZeroLagXcorr2Dtmem<ValueType>);
    }
    if (dimensionEM.compare("2d") == 0 && type.compare("viscoemem") == 0) {
        return ZeroLagXcorrPtr(new ZeroLagXcorr2Dviscoemem<ValueType>);
    }
    if (dimensionEM.compare("2d") == 0 && type.compare("viscotmem") == 0) {
        return ZeroLagXcorrPtr(new ZeroLagXcorr2Dviscotmem<ValueType>);
    }
    // 3D
//     if (dimensionEM.compare("3d") == 0 && type.compare("emem") == 0) {
//         return ZeroLagXcorrPtr(new ZeroLagXcorr3Demem<ValueType>);
//     }
//     if (dimensionEM.compare("3d") == 0 && type.compare("viscoemem") == 0) {
//         COMMON_THROWEXCEPTION("3Dviscoemem convolution is not implemented yet.");
//     }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return ZeroLagXcorrPtr();
}

template class KITGPI::ZeroLagXcorr::FactoryEM<double>;
template class KITGPI::ZeroLagXcorr::FactoryEM<float>;
