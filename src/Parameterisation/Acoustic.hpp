
#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include "Parameterisation.hpp"
#include <PartitionedInOut/PartitionedInOut.hpp>

namespace KITGPI
{

    //! \brief Parameterisation namespace
    namespace Parameterisation
    {

        //! Class for Parameterisation for acoustic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the acoustic finite-difference simulation.
         */
        template <typename ValueType>
        class Acoustic : public Parameterisation<ValueType>
        {
          public:
            //! Default constructor.
            Acoustic(){};

            //! Destructor, releases all allocated resources.
            ~Acoustic(){};

            explicit Acoustic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            explicit Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus_const, scai::lama::Scalar rho_const);
            explicit Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn);

            //! Copy Constructor.
            Acoustic(const Acoustic &rhs);

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus_const, scai::lama::Scalar rho_const);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn) override;

            void write(std::string filename, IndexType partitionedOut) const override;

            /* Getter methods for not requiered parameters */
            scai::lama::Vector const &getVelocityS() override;
            scai::lama::Vector const &getTauP() override;
            scai::lama::Vector const &getTauS() override;
            IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;

            /* Overloading Operators */
            KITGPI::Parameterisation::Acoustic<ValueType> operator*(scai::lama::Scalar rhs);
            KITGPI::Parameterisation::Acoustic<ValueType> &operator*=(scai::lama::Scalar const &rhs);
            KITGPI::Parameterisation::Acoustic<ValueType> operator+(KITGPI::Parameterisation::Acoustic<ValueType> const &rhs);
            KITGPI::Parameterisation::Acoustic<ValueType> &operator+=(KITGPI::Parameterisation::Acoustic<ValueType> const &rhs);
            KITGPI::Parameterisation::Acoustic<ValueType> operator-(KITGPI::Parameterisation::Acoustic<ValueType> const &rhs);
            KITGPI::Parameterisation::Acoustic<ValueType> &operator-=(KITGPI::Parameterisation::Acoustic<ValueType> const &rhs);
            KITGPI::Parameterisation::Acoustic<ValueType> &operator=(KITGPI::Parameterisation::Acoustic<ValueType> const &rhs);

            void minusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
            void plusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
            void assign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);

          private:
            using Parameterisation<ValueType>::density;
            using Parameterisation<ValueType>::velocityP;

            /* Not requiered parameters */
            using Parameterisation<ValueType>::velocityS;
            using Parameterisation<ValueType>::tauP;
            using Parameterisation<ValueType>::tauS;
            using Parameterisation<ValueType>::relaxationFrequency;
            using Parameterisation<ValueType>::numRelaxationMechanisms;
        };
    }
}
