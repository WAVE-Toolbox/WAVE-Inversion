#pragma once

#include <scai/lama.hpp>
#include <iostream>

#include <Acquisition/Receivers.hpp>
#include <Wavefields/Wavefields.hpp>
#include "../../Gradient/Gradient.hpp"

namespace KITGPI
{
    
    //! \brief Gradient preconditioning namespace
    namespace Preconditioning
    {
        /*! \brief Class to precondition a gradient based on migration weight K1 from Plessix and Mulder, 2004.
         * 
         *  It offers the possibility to calculate and apply the inverse of an approximation of the diagonal of the Hessian. 
         *  For this, the forward propagated wavefield (which is used in the misfit) is squared and integrated for all time steps.
         *  If applied to the gradient, source artifacts will be reduced and badly illuminated model parts will be enhanced. 
         *  In the current implementation, the preconditioning can only be applied to a single shot!
         */
        template <typename ValueType>
        class EnergyPreconditioning
        {

        public:
            
            /* Default constructor and destructor */
            EnergyPreconditioning(){};
            ~EnergyPreconditioning(){};

            void init(scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config);
            void intSquaredWavefields(KITGPI::Wavefields::Wavefields<ValueType> &wavefield, ValueType DT); //!< Integrate squared wavefields 
            void apply(KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, scai::IndexType shotNumber);
            void resetApproxHessian();
            
        private:
            scai::lama::DenseVector<ValueType> approxHessian;            // approximation of the diagonal of the inverse of the Hessian 
            scai::lama::DenseVector<ValueType> wavefieldSeismogramType;
//             scai::lama::DenseVector<ValueType> approxRecGreensFct;       // approximation of the receiver Green's function
            bool saveApproxHessian;
            std::string approxHessianName;
            ValueType epsilonHessian;

            
        };
    }
}
