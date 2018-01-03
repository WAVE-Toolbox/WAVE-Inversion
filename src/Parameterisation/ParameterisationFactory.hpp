

#pragma once

#include "Acoustic.hpp"
#include "Elastic.hpp"
#include "Parameterisation.hpp"
#include "Viscoelastic.hpp"
#include <string>

namespace KITGPI
{

    namespace Parameterisation
    {

        //! \brief Factory class.
        template <typename ValueType>
        class Factory
        {

          public:
            //! \brief Declare Parameterisation pointer
            typedef typename Parameterisation<ValueType>::ParameterisationPtr ParameterisationPtr;

            Factory() = delete;
            Factory(Factory const &) = delete;
            void operator=(Factory const &) = delete;

            /*! \brief Create the right simmulation with factory methode.
             *
             \param type Simmulation type (acoustic, elsstic, viscoelastic)
             */
            static ParameterisationPtr Create(std::string type);
        };
    }
}
