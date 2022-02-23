#pragma once

#include <scai/lama.hpp>
#include <iostream>

#include <IO/IO.hpp>
#include "../Gradient/Gradient.hpp"
#include <Acquisition/Receivers.hpp>
#include <Wavefields/Wavefields.hpp>

namespace KITGPI
{
    
    //! \brief Gradient preconditioning namespace
    namespace Preconditioning
    {
        /*! \brief Class to precondition a gradient based on migration weight K1 from Plessix and Mulder, 2004 (see also Shin et al., 2001) for useEnergyPreconditioning = 1. useEnergyPreconditioning = 2 is used for preconditioning source and receiver positions (see Kurzmann et al., 2013).
         * 
         *  It offers the possibility to calculate and apply the inverse of an approximation of the diagonal of the Hessian. 
         *  For this, the forward propagated velocity wavefield (all components available) is squared and integrated for all time steps.
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
            void intSquaredWavefields(KITGPI::Wavefields::Wavefields<ValueType> &wavefield, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldAdjoint, ValueType DT); //!< Integrate squared wavefields 
            
            void resetApproxHessian();
            
            void apply(KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, scai::IndexType shotNumber, scai::IndexType fileFormat);
            void applyTransform(scai::lama::Matrix<ValueType> const &lhs);
            
        private:
            scai::lama::DenseVector<ValueType> approxHessian;            // approximation of the diagonal of the inverse of the Hessian 
            scai::lama::DenseVector<ValueType> approxHessianAdjoint; 
            scai::lama::DenseVector<ValueType> wavefieldX;
            scai::lama::DenseVector<ValueType> wavefieldY;
            scai::lama::DenseVector<ValueType> wavefieldZ;
            scai::IndexType useEnergyPreconditioning;
            bool saveApproxHessian;
            std::string approxHessianName;
            std::string dimension; 
            std::string equationType;
            ValueType epsilonHessian;     
            bool isSeismic;
        };
    }
}
