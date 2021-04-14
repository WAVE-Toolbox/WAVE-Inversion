

#pragma once

#include "EMEM.hpp"
#include "ViscoEMEM.hpp"
#include "Gradient.hpp"
#include <string>

namespace KITGPI
{

    namespace Gradient
    {

        //! \brief Gradient factory class.
        template <typename ValueType>
        class FactoryEM
        {

          public:
            //! \brief Declare Gradient pointer
            typedef typename GradientEM<ValueType>::GradientPtr GradientPtr;

            FactoryEM() = delete;
            FactoryEM(FactoryEM const &) = delete;
            void operator=(FactoryEM const &) = delete;

            /*! \brief Create the right simulation with factory method.
             *
             \param type Simulation type (elastic, viscoelastic)
             */
            static GradientPtr Create(std::string type);
        };
    }
}
