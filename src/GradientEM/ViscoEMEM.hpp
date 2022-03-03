
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

        /*! \brief Class to store the gradients for viscoemem inversion
         *
         */
        template <typename ValueType>
        class ViscoEMEM : public GradientEM<ValueType>
        {
          public:
            //! Default constructor.
            ViscoEMEM() { equationType = "viscoemem"; };

            //! Destructor, releases all allocated resources.
            ~ViscoEMEM(){};

            explicit ViscoEMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            //! Copy Constructor.
            ViscoEMEM(const ViscoEMEM &rhs);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType electricConductivity_const, ValueType dielectricPermittivity_const, ValueType tauElectricConductivity_const, ValueType tauDielectricPermittivity_const);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void resetGradient();

            void write(std::string filename, scai::IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const override;

            std::string getEquationType() const;

            /* Getter methods for not required parameters */

            void estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
            void calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPrioriEM, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            
            void scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config) override;
            void applyEnergyPreconditioning(ValueType epsilonHessian, scai::IndexType saveApproxHessian, std::string filename, scai::IndexType fileFormat) override;
            void normalize();
            
            /* Overloading Operators */
            KITGPI::Gradient::ViscoEMEM<ValueType> operator*(ValueType rhs);
            KITGPI::Gradient::ViscoEMEM<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Gradient::ViscoEMEM<ValueType> operator+(KITGPI::Gradient::ViscoEMEM<ValueType> const &rhs);
            KITGPI::Gradient::ViscoEMEM<ValueType> &operator+=(KITGPI::Gradient::ViscoEMEM<ValueType> const &rhs);
            KITGPI::Gradient::ViscoEMEM<ValueType> operator-(KITGPI::Gradient::ViscoEMEM<ValueType> const &rhs);
            KITGPI::Gradient::ViscoEMEM<ValueType> &operator-=(KITGPI::Gradient::ViscoEMEM<ValueType> const &rhs);
            KITGPI::Gradient::ViscoEMEM<ValueType> &operator=(KITGPI::Gradient::ViscoEMEM<ValueType> const &rhs);

            void minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void assign(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs);
            void timesAssign(ValueType const &rhs);
            void timesAssign(scai::lama::Vector<ValueType> const &rhs);

            void sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot);        
            
            void sumGradientPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D cutCoordinate) override;            
            void smooth(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config) override;
            void applyMedianFilter(scai::dmemo::CommunicatorPtr commAll, KITGPI::Configuration::Configuration config) override;
            
          private:
            using Gradient<ValueType>::equationType;

            using Gradient<ValueType>::electricConductivity;
            using Gradient<ValueType>::dielectricPermittivity;
            using Gradient<ValueType>::tauElectricConductivity;
            using Gradient<ValueType>::tauDielectricPermittivity;
            using Gradient<ValueType>::relaxationFrequency;
            using Gradient<ValueType>::numRelaxationMechanisms;
            using Gradient<ValueType>::porosity;
            using Gradient<ValueType>::saturation;
            using Gradient<ValueType>::reflectivity;
            
            using Gradient<ValueType>::weightingVector;
            using Gradient<ValueType>::workflowInner;
            using Gradient<ValueType>::GaussianKernel;
            using Gradient<ValueType>::PX;
            using Gradient<ValueType>::PY;
                        
            /* Not required parameters */
        };
    }
}
