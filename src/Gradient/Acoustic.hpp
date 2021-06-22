
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

#include "../Workflow/Workflow.hpp"
#include <Modelparameter/Acoustic.hpp>

namespace KITGPI
{

    //! \brief Gradient namespace
    namespace Gradient
    {

        /*! \brief Class to store the gradients for acoustic inversion
         *
         */
        template <typename ValueType>
        class Acoustic : public Gradient<ValueType>
        {
          public:
            //! Default constructor.
            Acoustic() { equationType = "acoustic"; };

            //! Destructor, releases all allocated resources.
            ~Acoustic(){};

            explicit Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            //! Copy Constructor.
            Acoustic(const Acoustic &rhs);

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType rho_const);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            /*! \brief Set all wavefields to zero.
             */
            void resetGradient()
            {
                this->resetParameter(velocityP);
                this->resetParameter(density);
            };

            void write(std::string filename, scai::IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const override;

            std::string getEquationType() const;

            /* Getter methods for not requiered parameters */
            scai::lama::Vector<ValueType> const &getVelocityS() override;
            scai::lama::Vector<ValueType> const &getTauP() override;
            scai::lama::Vector<ValueType> const &getTauS() override;
            scai::IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;

            void estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) override; 
            void calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) override; 
            void calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
                        
            
            void calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPriori, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
            void applyMedianFilter(KITGPI::Configuration::Configuration config) override;
            
            void scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config) override;
            void normalize();

            /* Overloading Operators */
            KITGPI::Gradient::Acoustic<ValueType> operator*(ValueType rhs);
            KITGPI::Gradient::Acoustic<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Gradient::Acoustic<ValueType> operator+(KITGPI::Gradient::Acoustic<ValueType> const &rhs);
            KITGPI::Gradient::Acoustic<ValueType> &operator+=(KITGPI::Gradient::Acoustic<ValueType> const &rhs);
            KITGPI::Gradient::Acoustic<ValueType> operator-(KITGPI::Gradient::Acoustic<ValueType> const &rhs);
            KITGPI::Gradient::Acoustic<ValueType> &operator-=(KITGPI::Gradient::Acoustic<ValueType> const &rhs);
            KITGPI::Gradient::Acoustic<ValueType> &operator=(KITGPI::Gradient::Acoustic<ValueType> const &rhs);

            friend KITGPI::Modelparameter::Acoustic<ValueType> operator-(KITGPI::Modelparameter::Acoustic<ValueType> const &lhs, KITGPI::Gradient::Acoustic<ValueType> const &rhs)
            {
                KITGPI::Modelparameter::Acoustic<ValueType> result(lhs);
                result -= rhs;
                return result;
            };

            friend KITGPI::Modelparameter::Acoustic<ValueType> const &operator-=(KITGPI::Modelparameter::Acoustic<ValueType> &lhs, KITGPI::Gradient::Acoustic<ValueType> const &rhs)
            {
                scai::lama::DenseVector<ValueType> temp;
                temp = lhs.getVelocityP() - rhs.velocityP;
                lhs.setVelocityP(temp);
                temp = lhs.getDensity() - rhs.density;
                lhs.setDensity(temp);
                return lhs;
            };

            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void assign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void timesAssign(ValueType const &rhs);
            void timesAssign(scai::lama::Vector<ValueType> const &rhs);

            void sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot);
            
            void sumGradientPerShot(KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth) override;
            
          private:
            using Gradient<ValueType>::equationType;

            using Gradient<ValueType>::density;
            using Gradient<ValueType>::velocityP;
            using Gradient<ValueType>::porosity;
            using Gradient<ValueType>::saturation;
            using Gradient<ValueType>::workflowInner;

            /* Not requiered parameters */
            using Gradient<ValueType>::velocityS;
            using Gradient<ValueType>::tauP;
            using Gradient<ValueType>::tauS;
            using Gradient<ValueType>::relaxationFrequency;
            using Gradient<ValueType>::numRelaxationMechanisms;
        };
    }
}
