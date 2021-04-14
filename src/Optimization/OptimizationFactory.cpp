#include "OptimizationFactory.hpp"

template <typename ValueType>
typename KITGPI::Optimization::Optimization<ValueType>::OptimizationPtr KITGPI::Optimization::Factory<ValueType>::Create(std::string type)
{
    
    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(type.compare("steepestdescent") == 0 || type.compare("conjugategradient") == 0 || type.compare("lbfgs") == 0 || type.compare("truncatednewton") == 0, "Unknown type");

    if (type.compare("steepestdescent") == 0) {
        return OptimizationPtr(new SteepestDescent<ValueType>);
    }
    if (type.compare("conjugategradient") == 0) {
        return OptimizationPtr(new ConjugateGradient<ValueType>);
    }
    if (type.compare("lbfgs") == 0) {
        //     return OptimizationPtr(new LBFGS<ValueType>);
        COMMON_THROWEXCEPTION("No LBFGS implemented");
    }
    if (type.compare("truncatednewton") == 0) {
//         return OptimizationPtr(new SteepestDescent<ValueType>);
        COMMON_THROWEXCEPTION("No truncatednewton implemented");
    }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return OptimizationPtr();
}

template class KITGPI::Optimization::Factory<float>;
template class KITGPI::Optimization::Factory<double>;
