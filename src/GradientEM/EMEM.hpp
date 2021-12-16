
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

#include "GradientEM.hpp"

namespace KITGPI
{

    //! \brief Gradient namespace
    namespace Gradient
    {

        /*! \brief Class to store the gradients for emem inversion
         *
         */
        template <typename ValueType>
        class EMEM : public GradientEM<ValueType>
        {
          public:
            //! Default constructor.
            EMEM() { equationType = "emem"; };

            //! Destructor, releases all allocated resources.
            ~EMEM(){};

            explicit EMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            //! Copy Constructor.
            EMEM(const EMEM &rhs);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType electricConductivity_const, ValueType dielectricPermittivity_const);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void resetGradient();

            void write(std::string filename, scai::IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const override;

            std::string getEquationType() const;

            /* Getter methods for not required parameters */
            scai::lama::Vector<ValueType> const &getTauElectricConductivity() override;
            scai::lama::Vector<ValueType> const &getTauDielectricPermittivity() override;
            scai::IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;

            void estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            void applyMedianFilter(KITGPI::Configuration::Configuration config, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
            void calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPrioriEM, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
            void scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config) override;
            void applyEnergyPreconditioning(ValueType epsilonHessian, scai::IndexType saveApproxHessian, std::string filename, scai::IndexType fileFormat) override;
            void normalize();
            
            /* Overloading Operators */
            KITGPI::Gradient::EMEM<ValueType> operator*(ValueType rhs);
            KITGPI::Gradient::EMEM<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Gradient::EMEM<ValueType> operator+(KITGPI::Gradient::EMEM<ValueType> const &rhs);
            KITGPI::Gradient::EMEM<ValueType> &operator+=(KITGPI::Gradient::EMEM<ValueType> const &rhs);
            KITGPI::Gradient::EMEM<ValueType> operator-(KITGPI::Gradient::EMEM<ValueType> const &rhs);
            KITGPI::Gradient::EMEM<ValueType> &operator-=(KITGPI::Gradient::EMEM<ValueType> const &rhs);
            KITGPI::Gradient::EMEM<ValueType> &operator=(KITGPI::Gradient::EMEM<ValueType> const &rhs);

            void minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void assign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void timesAssign(ValueType const &rhs);
            void timesAssign(scai::lama::Vector<ValueType> const &rhs);

            void sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot);
            
            void sumGradientPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth) override;
            void smoothGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType FCmax) override;
            
          private:
            using Gradient<ValueType>::equationType;

            using Gradient<ValueType>::electricConductivity;
            using Gradient<ValueType>::dielectricPermittivity;
            using Gradient<ValueType>::porosity;
            using Gradient<ValueType>::saturation;
            using Gradient<ValueType>::reflectivity;
            using Gradient<ValueType>::workflowInner;
                        
            /* Not required parameters */
            using Gradient<ValueType>::tauElectricConductivity;
            using Gradient<ValueType>::tauDielectricPermittivity;
            using Gradient<ValueType>::relaxationFrequency;
            using Gradient<ValueType>::numRelaxationMechanisms;
            using Gradient<ValueType>::weightingVector;
        };
    }
}
