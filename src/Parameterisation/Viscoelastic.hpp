
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

        //! Class for Parameterisation for visco-elastic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the visco-elastic finite-difference simulation.
         */
        template <typename ValueType>
        class Viscoelastic : public Parameterisation<ValueType>
        {
          public:
            //! Default constructor.
            Viscoelastic(){};

            //! Destructor, releases all allocated resources.
            ~Viscoelastic(){};

            explicit Viscoelastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            explicit Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus_const, scai::lama::Scalar sWaveModulus_const, scai::lama::Scalar rho_const, scai::lama::Scalar tauP_const, scai::lama::Scalar tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            explicit Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn);

            //! Copy Constructor.
            Viscoelastic(const Viscoelastic &rhs);

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar velocityP_const, scai::lama::Scalar velocityS_const, scai::lama::Scalar rho_const, scai::lama::Scalar tauP_const, scai::lama::Scalar tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn) override;

            void initRelaxationMechanisms(IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);

            void write(std::string filename, IndexType partitionedOut) const override;

            void minusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
            void plusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
            void assign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);

            /* Overloading Operators */
            KITGPI::Parameterisation::Viscoelastic<ValueType> operator*(scai::lama::Scalar rhs);
            KITGPI::Parameterisation::Viscoelastic<ValueType> &operator*=(scai::lama::Scalar const &rhs);
            KITGPI::Parameterisation::Viscoelastic<ValueType> operator+(KITGPI::Parameterisation::Viscoelastic<ValueType> const &rhs);
            KITGPI::Parameterisation::Viscoelastic<ValueType> &operator+=(KITGPI::Parameterisation::Viscoelastic<ValueType> const &rhs);
            KITGPI::Parameterisation::Viscoelastic<ValueType> operator-(KITGPI::Parameterisation::Viscoelastic<ValueType> const &rhs);
            KITGPI::Parameterisation::Viscoelastic<ValueType> &operator-=(KITGPI::Parameterisation::Viscoelastic<ValueType> const &rhs);
            KITGPI::Parameterisation::Viscoelastic<ValueType> &operator=(KITGPI::Parameterisation::Viscoelastic<ValueType> const &rhs);

          private:
            using Parameterisation<ValueType>::density;
            using Parameterisation<ValueType>::velocityP;
            using Parameterisation<ValueType>::velocityS;

            using Parameterisation<ValueType>::tauS;
            using Parameterisation<ValueType>::tauP;
            using Parameterisation<ValueType>::relaxationFrequency;
            using Parameterisation<ValueType>::numRelaxationMechanisms;
        };
    }
}