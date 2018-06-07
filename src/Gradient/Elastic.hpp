
#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include "Gradient.hpp"
#include <PartitionedInOut/PartitionedInOut.hpp>
#include "../Workflow/Workflow.hpp"

namespace KITGPI
{

    //! \brief Gradient namespace
    namespace Gradient
    {

        /*! \brief Class to store the gradients for elastic inversion
         *
         */
        template <typename ValueType>
        class Elastic : public Gradient<ValueType>
        {
          public:
            //! Default constructor.
            Elastic(){};

            //! Destructor, releases all allocated resources.
            ~Elastic(){};

            explicit Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            //! Copy Constructor.
            Elastic(const Elastic &rhs);

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType pWaveModulus, ValueType sWaveModulus, ValueType rho);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;


            /*! \brief Set all wavefields to zero.
            */
            void resetGradient()
            {
                this->resetParameter(velocityP);
                this->resetParameter(velocityS);
                this->resetParameter(density);
            };

            void write(std::string filename, scai::IndexType partitionedOut, KITGPI::Workflow::Workflow<ValueType> const &workflow) const override;

            /* Getter methods for not requiered parameters */
            scai::lama::Vector<ValueType> const &getTauP() override;
            scai::lama::Vector<ValueType> const &getTauS() override;
            scai::IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;

            void estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow);

            /* Overloading Operators */
            KITGPI::Gradient::Elastic<ValueType> operator*(ValueType rhs);
            KITGPI::Gradient::Elastic<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Gradient::Elastic<ValueType> operator+(KITGPI::Gradient::Elastic<ValueType> const &rhs);
            KITGPI::Gradient::Elastic<ValueType> &operator+=(KITGPI::Gradient::Elastic<ValueType> const &rhs);
            KITGPI::Gradient::Elastic<ValueType> operator-(KITGPI::Gradient::Elastic<ValueType> const &rhs);
            KITGPI::Gradient::Elastic<ValueType> &operator-=(KITGPI::Gradient::Elastic<ValueType> const &rhs);
            KITGPI::Gradient::Elastic<ValueType> &operator=(KITGPI::Gradient::Elastic<ValueType> const &rhs);

            void minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void assign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void timesAssign(ValueType const &rhs);
            void timesAssign(scai::lama::Vector<ValueType> const &rhs);

          private:

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
