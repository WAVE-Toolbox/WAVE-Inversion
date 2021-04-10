#pragma once

#include <string>

#include "./Optimization.hpp"
#include "./SteepestDescent.hpp"
#include "./ConjugateGradient.hpp"

namespace KITGPI
{

    namespace Optimization
    {

        //! \brief Optimization factory class.
        template <typename ValueType>
        class Factory
        {

          public:
              
            //! \brief Declare Optimization pointer
            typedef typename Optimization<ValueType>::OptimizationPtr OptimizationPtr;

            Factory() = delete;
            Factory(Factory const &) = delete;
            void operator=(Factory const &) = delete;

            /*! \brief Create the right optimization with factory method.
             *
             \param type Optimization type (steepest descent, conjugate gradientEM, L-BFGS, etc.)
             */
            static OptimizationPtr Create(std::string type);
        };
    }
}
