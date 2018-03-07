#pragma once

#include <scai/lama.hpp>
#include <iostream>

#include <Acquisition/Receivers.hpp>
#include <Wavefields/Wavefields.hpp>
#include "../../Gradient/Gradient.hpp"

namespace KITGPI
{
    
    //! \brief Preconditioning namespace
    namespace Preconditioning
    {
        /*! \brief EnergyPreconditioning class (migration weight K1 after Plessix and Mulder, 2004)
         *  This class offers the possibility to calculate an approximation of the diagonal of the inverse of the Hessian.
         *  If applied to the gradient, source artifacts will be reduced and badly illuminated model parts will be enhanced. 
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
            void apply(KITGPI::Gradient::Gradient<ValueType> &gradient, IndexType shotNumber);
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
