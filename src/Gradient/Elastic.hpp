
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

#include "Gradient.hpp"
#include <PartitionedInOut/PartitionedInOut.hpp>

namespace KITGPI
{

    //! \brief Gradient namespace
    namespace Gradient
    {

        //! Class for Gradient for elastic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the elastic finite-difference simulation.
         */
        template <typename ValueType>
        class Elastic : public Gradient<ValueType>
        {
          public:
            //! Default constructor.
            Elastic(){};

            //! Destructor, releases all allocated resources.
            ~Elastic(){};

            explicit Elastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            //! Copy Constructor.
            Elastic(const Elastic &rhs);

            void init(Configuration::Configuration const &config,scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus, scai::lama::Scalar sWaveModulus, scai::lama::Scalar rho);
            void init(Configuration::Configuration const &config,scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;


            /*! \brief Set all wavefields to zero.
            */
            void resetGradient()
            {
                this->resetParameter(velocityP);
                this->resetParameter(velocityS);
                this->resetParameter(density);
            };

            void write(std::string filename, IndexType partitionedOut) const override;

            /* Getter methods for not requiered parameters */
            scai::lama::Vector const &getTauP() override;
            scai::lama::Vector const &getTauS() override;
            IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;

            void estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &/*correlatedWavefields*/, KITGPI::Modelparameter::Modelparameter<ValueType> const &/*model*/, ValueType /*DT*/) override
            {
                COMMON_THROWEXCEPTION("estimate is not implemented for viscoelastic gradients,yet ");
            };
            void scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model);

            /* Overloading Operators */
            KITGPI::Gradient::Elastic<ValueType> operator*(scai::lama::Scalar rhs);
            KITGPI::Gradient::Elastic<ValueType> &operator*=(scai::lama::Scalar const &rhs);
            KITGPI::Gradient::Elastic<ValueType> operator+(KITGPI::Gradient::Elastic<ValueType> const &rhs);
            KITGPI::Gradient::Elastic<ValueType> &operator+=(KITGPI::Gradient::Elastic<ValueType> const &rhs);
            KITGPI::Gradient::Elastic<ValueType> operator-(KITGPI::Gradient::Elastic<ValueType> const &rhs);
            KITGPI::Gradient::Elastic<ValueType> &operator-=(KITGPI::Gradient::Elastic<ValueType> const &rhs);
            KITGPI::Gradient::Elastic<ValueType> &operator=(KITGPI::Gradient::Elastic<ValueType> const &rhs);

            void minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void assign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void timesAssign(scai::lama::Scalar const &rhs);
            void timesAssign(scai::lama::Vector const &rhs);

          private:
            using Gradient<ValueType>::invertForVp;
            using Gradient<ValueType>::invertForVs;
            using Gradient<ValueType>::invertForDensity;

            using Gradient<ValueType>::density;
            using Gradient<ValueType>::velocityP;
            using Gradient<ValueType>::velocityS;

            /* Not requiered parameters */
            using Gradient<ValueType>::tauP;
            using Gradient<ValueType>::tauS;
            using Gradient<ValueType>::relaxationFrequency;
            using Gradient<ValueType>::numRelaxationMechanisms;
        };
    }
}
