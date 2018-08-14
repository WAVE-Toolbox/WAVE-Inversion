#pragma once

#include "ZeroLagXcorr.hpp"
#include "ZeroLagXcorr2Dacoustic.hpp"
#include "ZeroLagXcorr2Delastic.hpp"
//#include "ZeroLagXcorr2Dvisco.hpp"
#include "ZeroLagXcorr3Dacoustic.hpp"
#include "ZeroLagXcorr3Delastic.hpp"
//#include "ZeroLagXcorr3Dvisco.hpp"

namespace KITGPI
{

    //! \brief ZeroLagXcorr namespace
    namespace ZeroLagXcorr
    {

        //! \brief ZeroLagXcorr factory class.
        template <typename ValueType>
        class Factory
        {

          public:
            //! \brief Declare ZeroLagXcorr pointer
            typedef typename ZeroLagXcorr<ValueType>::ZeroLagXcorrPtr ZeroLagXcorrPtr;

            Factory() = delete;
            Factory(Factory const &) = delete;
            void operator=(Factory const &) = delete;

            /*! \brief Create the right simmulation with factory methode.
             *
             \param dimension Dimension of the model (2D, 3D)
             \param type Simmulation type (acoustic, elsstic, viscoelastic)
             */
            static ZeroLagXcorrPtr Create(std::string dimension, std::string type);
        };
    }
}
