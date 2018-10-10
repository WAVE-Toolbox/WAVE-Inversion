
#pragma once

#include "Taper.hpp"
#include "Taper2D.hpp"
#include "Taper1D.hpp"


namespace KITGPI
{

    namespace Taper
    {

        //! \brief Factory class.
        template <typename ValueType>
        class Factory
        {
          public:
            //! \brief Declare Taper pointer
            typedef typename Taper<ValueType>::TaperPtr TaperPtr;

            Factory() = delete;
            Factory(Factory const &) = delete;
            void operator=(Factory const &) = delete;

            /*! \brief Create the right taper with factory methode.
             *
             \param dimension Dimension of the taper (1D, 2D)
             \param type Taper type (linear, cosine)
             */
            static TaperPtr Create(std::string dimension);
        };
    } /* end namespace Taper */
} /* end namespace KITGPI */