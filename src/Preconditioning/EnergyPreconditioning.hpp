#pragma once

#include <scai/lama.hpp>
#include <iostream>

#include <IO/IO.hpp>
#include <Acquisition/Receivers.hpp>
#include <Wavefields/Wavefields.hpp>
#include "../Gradient/Gradient.hpp"
#include <AcquisitionEM/Receivers.hpp>
#include <WavefieldsEM/Wavefields.hpp>
#include "../GradientEM/Gradient.hpp"

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

            void init(scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config, bool isSeismic);
            void intSquaredWavefields(KITGPI::Wavefields::Wavefields<ValueType> &wavefield, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldBack, ValueType DT); //!< Integrate squared wavefields 
            void apply(KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, scai::IndexType shotNumber, scai::IndexType fileFormat);
            
            void intSquaredWavefields(KITGPI::Wavefields::WavefieldsEM<ValueType> &wavefield, KITGPI::Wavefields::WavefieldsEM<ValueType> &wavefieldBack, ValueType DT); //!< Integrate squared wavefieldsEM 
            void apply(KITGPI::Gradient::GradientEM<ValueType> &gradientPerShotEM, scai::IndexType shotNumber, scai::IndexType fileFormat);
            void resetApproxHessian();
            
        private:
            scai::lama::DenseVector<ValueType> approxHessian;            // approximation of the diagonal of the inverse of the Hessian 
            scai::lama::DenseVector<ValueType> wavefieldVX;
            scai::lama::DenseVector<ValueType> wavefieldVY;
            scai::lama::DenseVector<ValueType> wavefieldVZ;
            scai::lama::DenseVector<ValueType> wavefieldEX;
            scai::lama::DenseVector<ValueType> wavefieldEY;
            scai::lama::DenseVector<ValueType> wavefieldEZ;
            scai::IndexType useEnergyPreconditioning;
            bool saveApproxHessian;
            std::string approxHessianName;
            std::string dimension; 
            std::string equationType;
            ValueType epsilonHessian;            
        };
    }
}
