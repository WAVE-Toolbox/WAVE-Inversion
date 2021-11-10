#pragma once

#include "../Gradient/Gradient.hpp"

namespace KITGPI
{

    //! \brief Gradient namespace
    namespace Gradient
    {
        
        /*! \brief Abstract class to store gradients for inversion
         *
	     * This class implements some methods which are required by all derived classes.
         * As this class is an abstract class, all constructors are protected.
         */
        template <typename ValueType>
        class GradientEM : public Gradient<ValueType>
        {
          public:
            //! Default constructor.
            GradientEM(){};

            //! Default destructor.
            ~GradientEM(){};

            /* Common */
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            virtual void resetGradient() = 0;

            virtual void estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;
            
            virtual void calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) override; 
            virtual void calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) override; 
            virtual void calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) override;
            virtual ValueType calcCrossGradientMisfit() override;
            
            virtual void applyMedianFilter(KITGPI::Configuration::Configuration config, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

            virtual void calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPriori, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;
            
            virtual void write(std::string filename, scai::IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const = 0;

            virtual std::string getEquationType() const = 0;

            virtual void scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config) = 0;
            virtual void applyEnergyPreconditioning(ValueType epsilonHessian, scai::IndexType saveApproxHessian, std::string filename, scai::IndexType fileFormat) = 0;
            virtual void normalize() = 0;

            virtual void minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs) = 0;
            virtual void plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs) = 0;
            virtual void assign(KITGPI::Gradient::Gradient<ValueType> const &rhs) = 0;
            virtual void timesAssign(ValueType const &rhs) = 0;
            virtual void timesAssign(scai::lama::Vector<ValueType> const &rhs) = 0;

            virtual void sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot) = 0;
            
            virtual void sumGradientPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth) = 0;
            
            virtual void setInvertForParameters(std::vector<bool> setInvertForParameters) override;
            
            virtual void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs) = 0;

            /* Seismic */
            virtual scai::lama::DenseVector<ValueType> getBiotCoefficientDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) override; 
            virtual scai::lama::DenseVector<ValueType> getDensityDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) override;
            virtual scai::lama::DenseVector<ValueType> getMu_satDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) override;
            virtual scai::lama::DenseVector<ValueType> getK_satDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) override;
            virtual scai::lama::DenseVector<ValueType> getDensityDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) override;
            virtual scai::lama::DenseVector<ValueType> getK_satDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) override;
            
            virtual scai::lama::Vector<ValueType> const &getDensity() override;
            virtual scai::lama::Vector<ValueType> const &getDensity() const override;
            virtual scai::lama::Vector<ValueType> const &getVelocityP() override;
            virtual scai::lama::Vector<ValueType> const &getVelocityP() const override;
            virtual scai::lama::Vector<ValueType> const &getVelocityS() override;
            virtual scai::lama::Vector<ValueType> const &getVelocityS() const override;

            virtual scai::lama::Vector<ValueType> const &getTauP() override;
            virtual scai::lama::Vector<ValueType> const &getTauP() const override;
            virtual scai::lama::Vector<ValueType> const &getTauS() override;
            virtual scai::lama::Vector<ValueType> const &getTauS() const override;

            virtual void setDensity(scai::lama::Vector<ValueType> const &setDensity) override;
            virtual void setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP) override;
            virtual void setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS) override;

            virtual void setTauP(scai::lama::Vector<ValueType> const &setTauP) override;
            virtual void setTauS(scai::lama::Vector<ValueType> const &setTauS) override;

            /* EM */
            virtual scai::lama::Vector<ValueType> const &getElectricConductivity() override;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivity() const override;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivity() override;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivity() const override;
            
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivity() override;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivity() const override;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivity() override;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivity() const override;
            
            virtual void setElectricConductivity(scai::lama::Vector<ValueType> const &setElectricConductivity) override;
            virtual void setDielectricPermittivity(scai::lama::Vector<ValueType> const &setDielectricPermittivity) override; 
            
            virtual void setTauElectricConductivity(scai::lama::Vector<ValueType> const &setTauElectricConductivity) override;
            virtual void setTauDielectricPermittivity(scai::lama::Vector<ValueType> const &setTauDielectricPermittivity) override;
            
            virtual void applyParameterisation(ValueType &modelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation) override;
            virtual void applyParameterisation(scai::lama::DenseVector<ValueType> &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation) override;
            virtual void deleteParameterisation(scai::lama::DenseVector<ValueType> &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation) override;
            virtual void gradientParameterisation(scai::lama::DenseVector<ValueType> &vecGradientParameter, scai::lama::DenseVector<ValueType> vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation) override;
            virtual scai::lama::DenseVector<ValueType> getElectricConductivityDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) override;
            virtual scai::lama::DenseVector<ValueType> getDielectricPermittiviyDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) override;
            virtual scai::lama::DenseVector<ValueType> getElectricConductivityDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) override;
            virtual scai::lama::DenseVector<ValueType> getDielectricPermittiviyDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) override;

          protected:
            /* Common */
            using Gradient<ValueType>::equationType;

            using Gradient<ValueType>::numRelaxationMechanisms; //!< Number of relaxation mechanisms
            using Gradient<ValueType>::relaxationFrequency;           //!< Relaxation Frequency

            using Gradient<ValueType>::porosity; //!< Vector storing porosity.
            using Gradient<ValueType>::saturation; //!< Vector storing saturation
            using Gradient<ValueType>::reflectivity; //!< Vector storing reflectivity

            using Gradient<ValueType>::normalizeGradient;
            using Gradient<ValueType>::workflowInner;
            using Gradient<ValueType>::weightingVector;
            
            /* Seismic */
            using Gradient<ValueType>::density; //!< Vector storing Density.

            using Gradient<ValueType>::velocityP; //!< Vector storing P-wave velocity.
            using Gradient<ValueType>::velocityS; //!< Vector storing S-wave velocity.

            using Gradient<ValueType>::tauP; //!< Vector storing tauP for visco-elastic modelling.
            using Gradient<ValueType>::tauS; //!< Vector storing tauS for visco-elastic modelling.

            /* EM */
            using Gradient<ValueType>::electricConductivity; //!< Vector storing electricConductivity.
            using Gradient<ValueType>::dielectricPermittivity; //!< Vector storing dielectricPermittivity.
            
            using Gradient<ValueType>::tauElectricConductivity; //!< Vector storing tauElectricConductivity for visco-emem modelling.
            using Gradient<ValueType>::tauDielectricPermittivity; //!< Vector storing TauDielectricPermittivity for visco-emem modelling.    
        };
    }
}
