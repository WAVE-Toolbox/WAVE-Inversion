#pragma once

#include "Acoustic.hpp"
#include "Elastic.hpp"
#include "SH.hpp"
#include "ViscoSH.hpp"
#include "Viscoelastic.hpp"

#include "../GradientEM/EMEM.hpp"
#include "../GradientEM/ViscoEMEM.hpp"

#include <string>

namespace KITGPI
{

    namespace Gradient
    {

        //! \brief Gradient factory class.
        template <typename ValueType>
        class Factory
        {

          public:
            //! \brief Declare Gradient pointer
            typedef typename Gradient<ValueType>::GradientPtr GradientPtr;

            Factory() = delete;
            Factory(Factory const &) = delete;
            void operator=(Factory const &) = delete;

            /*! \brief Create the right simmulation with factory methode.
             *
             \param type Simmulation type (acoustic, elsstic, viscoelastic)
             */
            static GradientPtr Create(std::string type);
        };
    }
}
