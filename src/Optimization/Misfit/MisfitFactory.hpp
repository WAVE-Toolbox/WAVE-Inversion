#pragma once

#include <string>
#include "Misfit.hpp"
#include "MisfitL2.hpp"

namespace KITGPI
{

    namespace Misfit
    {

        //! \brief Misfit factory class.
        template <typename ValueType>
        class Factory
        {

          public:
            //! \brief Declare misfit pointer
            typedef typename Misfit<ValueType>::MisfitPtr MisfitPtr;

            Factory() = delete;
            Factory(Factory const &) = delete;
            void operator=(Factory const &) = delete;

            /*! \brief Create the right misfit with factory method.
             *
             \param type Misfit type (L2, etc.)
             */
            static MisfitPtr Create(std::string type);
        };
    }
}
