
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

#include "../Workflow/Workflow.hpp"
#include "Gradient.hpp"

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
            Elastic() { equationType = "elastic"; };

            //! Destructor, releases all allocated resources.
            ~Elastic(){};

            explicit Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            //! Copy Constructor.
            Elastic(const Elastic &rhs);

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho_const, ValueType porosity_const, ValueType saturation_const);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            /*! \brief Set all wavefields to zero.
            */
            void resetGradient()
            {
                this->resetParameter(velocityP);
                this->resetParameter(velocityS);
                this->resetParameter(density);
                this->resetParameter(porosity);
                this->resetParameter(saturation);
            };

            void write(std::string filename, scai::IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const override;

            std::string getEquationType() const;

            /* Getter methods for not requiered parameters */
            scai::lama::Vector<ValueType> const &getTauP() override;
            scai::lama::Vector<ValueType> const &getTauS() override;
            scai::IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;

            void estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> gradientTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) override; 
            void calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> gradientTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) override; 
            void calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> gradientTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
            
            void calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPriori, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
            void applyMedianFilter(KITGPI::Configuration::Configuration config) override;
            
            void scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config) override;
            void normalize();

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

            void sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot);
            
            void sumGradientPerShot(KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth) override;

          private:
            using Gradient<ValueType>::equationType;

            using Gradient<ValueType>::density;
            using Gradient<ValueType>::velocityP;
            using Gradient<ValueType>::velocityS;
            using Gradient<ValueType>::porosity;
            using Gradient<ValueType>::saturation;

            /* Not requiered parameters */
            using Gradient<ValueType>::tauP;
            using Gradient<ValueType>::tauS;
            using Gradient<ValueType>::relaxationFrequency;
            using Gradient<ValueType>::numRelaxationMechanisms;
        };
    }
}
