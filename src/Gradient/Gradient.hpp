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

#include "../Common/HostPrint.hpp"
#include <Configuration/Configuration.hpp>

#include "../Workflow/Workflow.hpp"
#include "../ZeroLagCrossCorrelation/ZeroLagXcorr.hpp"
#include <Common/Common.hpp>
#include <Modelparameter/Modelparameter.hpp>

#include <stdlib.h>
#include <iomanip>
#include "../Common/Common.hpp"
#include "../Misfit/Misfit.hpp"
#include <IO/IO.hpp>
#include "../Taper/Taper2D.hpp"

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
        class Gradient
        {
          public:
            //! Default constructor.
            Gradient() : numRelaxationMechanisms(0), normalizeGradient(false){};

            //! Default destructor.
            ~Gradient(){};

            //! \brief Gradient pointer
            typedef std::shared_ptr<Gradient<ValueType>> GradientPtr;

            /* Common */
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            virtual void resetGradient() = 0;

            virtual void estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;
            virtual void calcModelDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0; 
            virtual void calcCrossGradient(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0; 
            virtual void calcCrossGradientDerivative(KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivativesEM, KITGPI::Configuration::Configuration configEM, KITGPI::Taper::Taper2D<ValueType> modelTaper2DJoint, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;
            virtual ValueType calcCrossGradientMisfit() = 0;            
            
            virtual void applyMedianFilter(KITGPI::Configuration::Configuration config, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

            scai::lama::DenseVector<ValueType> calcStabilizingFunctionalGradientPerModel(scai::lama::DenseVector<ValueType> modelResidualVec, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit);
            virtual void calcStabilizingFunctionalGradient(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelPriori, KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;
            
            virtual void write(std::string filename, scai::IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const = 0;

            virtual std::string getEquationType() const = 0;

            virtual scai::IndexType getNumRelaxationMechanisms() const;
            virtual ValueType getRelaxationFrequency() const;
            bool getNormalizeGradient() const;

            virtual void setNumRelaxationMechanisms(scai::IndexType const setNumRelaxationMechanisms);
            virtual void setRelaxationFrequency(ValueType const setRelaxationFrequency);

            virtual scai::lama::Vector<ValueType> const &getPorosity();
            virtual scai::lama::Vector<ValueType> const &getPorosity() const;
            virtual scai::lama::Vector<ValueType> const &getSaturation();
            virtual scai::lama::Vector<ValueType> const &getSaturation() const;
            virtual scai::lama::Vector<ValueType> const &getReflectivity();
            virtual scai::lama::Vector<ValueType> const &getReflectivity() const;

            virtual void setPorosity(scai::lama::Vector<ValueType> const &setPorosity);
            virtual void setSaturation(scai::lama::Vector<ValueType> const &setSaturation);
            virtual void setReflectivity(scai::lama::Vector<ValueType> const &setReflectivity); 

            void setNormalizeGradient(bool const &normGrad);

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
            
            void calcWeightingVector(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, std::vector<scai::IndexType> uniqueShotInds);
            scai::lama::DenseVector<ValueType> getWeightingVector();

            virtual void setInvertForParameters(std::vector<bool> setInvertForParameters) = 0;
            std::vector<bool> getInvertForParameters();
            
            /* Operator overloading */
            /*lhs Base rhs Base */
            KITGPI::Gradient::Gradient<ValueType> &operator=(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            KITGPI::Gradient::Gradient<ValueType> &operator-=(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            KITGPI::Gradient::Gradient<ValueType> &operator+=(KITGPI::Gradient::Gradient<ValueType> const &rhs);
            KITGPI::Gradient::Gradient<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Gradient::Gradient<ValueType> &operator*=(scai::lama::Vector<ValueType> const &rhs);

            /*lhs: fd-Model-Base rhs: gradient Base */
            friend KITGPI::Modelparameter::Modelparameter<ValueType> &operator-=(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> &rhs)
            {
                rhs.minusAssign(lhs, rhs);
                return lhs;
            };

            virtual void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs) = 0;

            /* Seismic */
            virtual scai::lama::DenseVector<ValueType> getBiotCoefficientDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0; 
            virtual scai::lama::DenseVector<ValueType> getDensityDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0;
            virtual scai::lama::DenseVector<ValueType> getMu_satDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0;
            virtual scai::lama::DenseVector<ValueType> getK_satDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0;
            virtual scai::lama::DenseVector<ValueType> getDensityDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0;
            virtual scai::lama::DenseVector<ValueType> getK_satDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0;
            
            virtual scai::lama::Vector<ValueType> const &getDensity() = 0;
            virtual scai::lama::Vector<ValueType> const &getDensity() const = 0;
            virtual scai::lama::Vector<ValueType> const &getVelocityP() = 0;
            virtual scai::lama::Vector<ValueType> const &getVelocityP() const = 0;
            virtual scai::lama::Vector<ValueType> const &getVelocityS() = 0;
            virtual scai::lama::Vector<ValueType> const &getVelocityS() const = 0;

            virtual scai::lama::Vector<ValueType> const &getTauP() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauP() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauS() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauS() const = 0;

            virtual void setDensity(scai::lama::Vector<ValueType> const &setDensity) = 0;
            virtual void setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP) = 0;
            virtual void setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS) = 0;

            virtual void setTauP(scai::lama::Vector<ValueType> const &setTauP) = 0;
            virtual void setTauS(scai::lama::Vector<ValueType> const &setTauS) = 0;

            /* EM */
            virtual scai::lama::Vector<ValueType> const &getElectricConductivity() = 0;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivity() const = 0;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivity() = 0;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivity() const = 0;
            
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivity() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivity() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivity() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivity() const = 0;
            
            virtual void setElectricConductivity(scai::lama::Vector<ValueType> const &setElectricConductivity) = 0;
            virtual void setDielectricPermittivity(scai::lama::Vector<ValueType> const &setDielectricPermittivity) = 0; 
            
            virtual void setTauElectricConductivity(scai::lama::Vector<ValueType> const &setTauElectricConductivity) = 0;
            virtual void setTauDielectricPermittivity(scai::lama::Vector<ValueType> const &setTauDielectricPermittivity) = 0;
            
            virtual void applyParameterisation(ValueType &modelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation) = 0;
            virtual void applyParameterisation(scai::lama::DenseVector<ValueType> &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation) = 0;
            virtual void deleteParameterisation(scai::lama::DenseVector<ValueType> &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation) = 0;
            virtual void gradientParameterisation(scai::lama::DenseVector<ValueType> &vecGradientParameter, scai::lama::Vector<ValueType> const &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation) = 0;
            virtual scai::lama::DenseVector<ValueType> getElectricConductivityDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0;
            virtual scai::lama::DenseVector<ValueType> getDielectricPermittiviyDePorosity(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0;
            virtual scai::lama::DenseVector<ValueType> getElectricConductivityDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0;
            virtual scai::lama::DenseVector<ValueType> getDielectricPermittiviyDeSaturation(KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0;
            
          protected:
            /* Common */
            void initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value);
            void initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            void writeParameterisation(scai::lama::Vector<ValueType> const &vector, std::string filename, scai::IndexType fileFormat) const;

            void resetParameter(scai::lama::DenseVector<ValueType> &vector) { vector.setScalar(0); }

            std::string equationType;

            scai::IndexType numRelaxationMechanisms; //!< Number of relaxation mechanisms
            ValueType relaxationFrequency;           //!< Relaxation Frequency

            scai::lama::DenseVector<ValueType> porosity; //!< Vector storing porosity.
            scai::lama::DenseVector<ValueType> saturation; //!< Vector storing saturation
            scai::lama::DenseVector<ValueType> reflectivity; //!< Vector storing reflectivity

            bool normalizeGradient;
            KITGPI::Workflow::Workflow<ValueType> workflowInner;
            scai::lama::DenseVector<ValueType> weightingVector;
            
            /* Seismic */
            scai::lama::DenseVector<ValueType> density; //!< Vector storing Density.

            scai::lama::DenseVector<ValueType> velocityP; //!< Vector storing P-wave velocity.
            scai::lama::DenseVector<ValueType> velocityS; //!< Vector storing S-wave velocity.

            scai::lama::DenseVector<ValueType> tauP; //!< Vector storing tauP for visco-elastic modelling.
            scai::lama::DenseVector<ValueType> tauS; //!< Vector storing tauS for visco-elastic modelling.

            /* EM */
            scai::lama::DenseVector<ValueType> electricConductivity; //!< Vector storing electricConductivity.
            scai::lama::DenseVector<ValueType> dielectricPermittivity; //!< Vector storing dielectricPermittivity.
            
            scai::lama::DenseVector<ValueType> tauElectricConductivity; //!< Vector storing tauElectricConductivity for visco-emem modelling.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivity; //!< Vector storing TauDielectricPermittivity for visco-emem modelling.
            
          private:
            void allocateParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void readParameterisation(scai::lama::Vector<ValueType> &vector, std::string filename, scai::IndexType fileFormat);
        };
    }
}
