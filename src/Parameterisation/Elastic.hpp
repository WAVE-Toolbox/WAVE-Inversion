
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

        //! Class for Parameterisation for elastic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the elastic finite-difference simulation.
         */
        template <typename ValueType>
        class Elastic : public Parameterisation<ValueType>
        {
          public:
            //! Default constructor.
            Elastic(){};

            //! Destructor, releases all allocated resources.
            ~Elastic(){};

            explicit Elastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            explicit Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus_const, scai::lama::Scalar sWaveModulus_const, scai::lama::Scalar rho);
            explicit Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn);

            //! Copy Constructor.
            Elastic(const Elastic &rhs);

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus, scai::lama::Scalar sWaveModulus, scai::lama::Scalar rho);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn) override;
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void write(std::string filename, IndexType partitionedOut) const override;

            /* Getter methods for not requiered parameters */
            scai::lama::Vector const &getTauP() override;
            scai::lama::Vector const &getTauS() override;
            IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;

            /* Overloading Operators */
            KITGPI::Parameterisation::Elastic<ValueType> operator*(scai::lama::Scalar rhs);
            KITGPI::Parameterisation::Elastic<ValueType> &operator*=(scai::lama::Scalar const &rhs);
            KITGPI::Parameterisation::Elastic<ValueType> operator+(KITGPI::Parameterisation::Elastic<ValueType> const &rhs);
            KITGPI::Parameterisation::Elastic<ValueType> &operator+=(KITGPI::Parameterisation::Elastic<ValueType> const &rhs);
            KITGPI::Parameterisation::Elastic<ValueType> operator-(KITGPI::Parameterisation::Elastic<ValueType> const &rhs);
            KITGPI::Parameterisation::Elastic<ValueType> &operator-=(KITGPI::Parameterisation::Elastic<ValueType> const &rhs);
            KITGPI::Parameterisation::Elastic<ValueType> &operator=(KITGPI::Parameterisation::Elastic<ValueType> const &rhs);

            void minusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
            void plusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
            void assign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);

          private:
            using Parameterisation<ValueType>::density;
            using Parameterisation<ValueType>::velocityP;
            using Parameterisation<ValueType>::velocityS;

            /* Not requiered parameters */
            using Parameterisation<ValueType>::tauP;
            using Parameterisation<ValueType>::tauS;
            using Parameterisation<ValueType>::relaxationFrequency;
            using Parameterisation<ValueType>::numRelaxationMechanisms;
        };
    }
}
