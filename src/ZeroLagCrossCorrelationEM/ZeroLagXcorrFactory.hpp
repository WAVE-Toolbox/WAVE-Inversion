#pragma once

#include "ZeroLagXcorr.hpp"
#include "ZeroLagXcorr2Demem.hpp"
#include "ZeroLagXcorr2Dtmem.hpp"
#include "ZeroLagXcorr2Dviscoemem.hpp"
#include "ZeroLagXcorr2Dviscotmem.hpp"

namespace KITGPI
{

    //! \brief ZeroLagXcorr namespace
    namespace ZeroLagXcorr
    {

        //! \brief ZeroLagXcorr factory class.
        template <typename ValueType>
        class FactoryEM
        {

          public:
            //! \brief Declare ZeroLagXcorr pointer
            typedef typename ZeroLagXcorrEM<ValueType>::ZeroLagXcorrPtr ZeroLagXcorrPtr;

            FactoryEM() = delete;
            FactoryEM(FactoryEM const &) = delete;
            void operator=(FactoryEM const &) = delete;

            /*! \brief Create the right simulation with factory method.
             *
             \param dimensionEM Dimension of the modelEM (2D, 3D)
             \param type Simulation type ( elastic, viscoelastic)
             */
            static ZeroLagXcorrPtr Create(std::string dimensionEM, std::string type);
        };
    }
}
