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

#include <AcquisitionEM/Sources.hpp>
#include <AcquisitionEM/Receivers.hpp>
#include <WavefieldsEM/WavefieldsFactory.hpp>
#include "../WorkflowEM/Workflow.hpp"
#include "../ZeroLagCrossCorrelationEM/ZeroLagXcorr.hpp"
#include <Common/Common.hpp>
#include <Modelparameter/Modelparameter.hpp>
#include <ModelparameterEM/Modelparameter.hpp>
#include <ForwardSolver/Derivatives/Derivatives.hpp>

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
        class GradientEM
        {
          public:
            //! Default constructor.
            GradientEM() : numRelaxationMechanisms(0), normalizeGradient(false){};

            //! Default destructor.
            ~GradientEM(){};

            //! \brief Gradient pointer
            typedef std::shared_ptr<GradientEM<ValueType>> GradientPtr;

            /*! \brief Abstract initialisation function
             * Standard initialisation function
             \param ctx Context
             \param distEM Distribution
             */
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM) = 0;

            virtual void resetGradient() = 0;

            virtual void estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, ValueType DT, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) = 0;
            virtual void applyMedianFilter(KITGPI::Configuration::Configuration configEM) = 0;
            
            void exParameterisation(ValueType &modelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation);
            void exParameterisation(scai::lama::DenseVector<ValueType> &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation);
            void deParameterisation(scai::lama::DenseVector<ValueType> &vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation);
            void gradientParameterisation(scai::lama::DenseVector<ValueType> &vecGradientParameter, scai::lama::DenseVector<ValueType> vecModelParameter, ValueType const modelParameterReference, scai::IndexType parameterisation);
            scai::lama::DenseVector<ValueType> getConductivityDePorosity(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM);
            scai::lama::DenseVector<ValueType> getDielectricPermittiviyDePorosity(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM);
            scai::lama::DenseVector<ValueType> getConductivityDeSaturation(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM);
            scai::lama::DenseVector<ValueType> getDielectricPermittiviyDeSaturation(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM);
            
            scai::lama::DenseVector<ValueType> calcStabilizingFunctionalGradientPerModel(scai::lama::DenseVector<ValueType> modelResidualVec, KITGPI::Configuration::Configuration configEM, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM);
            virtual void calcStabilizingFunctionalGradient(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelPrioriEM, KITGPI::Configuration::Configuration configEM, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) = 0;
                           
            /*! \brief Abstract write function
             *
             * Standard write function
             *
             \param filename filename to write parameters (endings will be added by derived classes)
             \param fileFormat format of output file
             */
            virtual void write(std::string filename, scai::IndexType fileFormat, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM) const = 0;

            virtual std::string getEquationType() const = 0;

            virtual scai::lama::Vector<ValueType> const &getConductivityEM();
            virtual scai::lama::Vector<ValueType> const &getConductivityEM() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEM();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEM() const;
            
            virtual scai::lama::Vector<ValueType> const &getTauConductivityEM();
            virtual scai::lama::Vector<ValueType> const &getTauConductivityEM() const;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityEM();
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityEM() const;
            virtual scai::IndexType getNumRelaxationMechanisms() const;
            virtual ValueType getRelaxationFrequency() const;
            
            virtual scai::lama::Vector<ValueType> const &getPorosity();
            virtual scai::lama::Vector<ValueType> const &getPorosity() const;
            virtual scai::lama::Vector<ValueType> const &getSaturation();
            virtual scai::lama::Vector<ValueType> const &getSaturation() const;

            bool getNormalizeGradient() const;

            virtual void setConductivityEM(scai::lama::Vector<ValueType> const &setConductivityEM);
            virtual void setDielectricPermittivityEM(scai::lama::Vector<ValueType> const &setDielectricPermittivityEM); 
            
            virtual void setTauConductivityEM(scai::lama::Vector<ValueType> const &setTauConductivityEM);
            virtual void setTauDielectricPermittivityEM(scai::lama::Vector<ValueType> const &setTauDielectricPermittivityEM);
            virtual void setNumRelaxationMechanisms(scai::IndexType const setNumRelaxationMechanisms);
            virtual void setRelaxationFrequency(ValueType const setRelaxationFrequency);

            virtual void setPorosity(scai::lama::Vector<ValueType> const &setPorosity);
            virtual void setSaturation(scai::lama::Vector<ValueType> const &setSaturation); 

            void setNormalizeGradient(bool const &normGrad);

            virtual void scale(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM, KITGPI::Configuration::Configuration configEM) = 0;
            virtual void normalize() = 0;

            virtual void minusAssign(KITGPI::Gradient::GradientEM<ValueType> const &rhs) = 0;
            virtual void plusAssign(KITGPI::Gradient::GradientEM<ValueType> const &rhs) = 0;
            virtual void assign(KITGPI::Gradient::GradientEM<ValueType> const &rhs) = 0;
            virtual void timesAssign(ValueType const &rhs) = 0;
            virtual void timesAssign(scai::lama::Vector<ValueType> const &rhs) = 0;

            virtual void sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot) = 0;
                        
            virtual void sumGradientPerShot(KITGPI::Gradient::GradientEM<ValueType> &gradientPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType shotInd, scai::IndexType boundaryWidth) = 0;
            
            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Declare Sparse-Matrix
            SparseFormat getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate);
            
            scai::lama::SparseVector<ValueType> getEraseVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth);

            /* Operator overloading */
            /*lhs Base rhs Base */
            KITGPI::Gradient::GradientEM<ValueType> &operator=(KITGPI::Gradient::GradientEM<ValueType> const &rhs);
            KITGPI::Gradient::GradientEM<ValueType> &operator-=(KITGPI::Gradient::GradientEM<ValueType> const &rhs);
            KITGPI::Gradient::GradientEM<ValueType> &operator+=(KITGPI::Gradient::GradientEM<ValueType> const &rhs);
            KITGPI::Gradient::GradientEM<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Gradient::GradientEM<ValueType> &operator*=(scai::lama::Vector<ValueType> const &rhs);

            /*lhs: fd-Model-Base rhs: gradientEM Base */
            friend KITGPI::Modelparameter::ModelparameterEM<ValueType> &operator-=(KITGPI::Modelparameter::ModelparameterEM<ValueType> &lhs, KITGPI::Gradient::GradientEM<ValueType> &rhs)
            {
                rhs.minusAssign(lhs, rhs);
                return lhs;
            };

            virtual void minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> &lhs, KITGPI::Gradient::GradientEM<ValueType> const &rhs) = 0;            

          protected:
            void resetParameter(scai::lama::DenseVector<ValueType> &vector) { vector.setScalar(0); }

            std::string equationTypeEM;

            scai::lama::DenseVector<ValueType> conductivityEM; //!< Vector storing conductivityEM.
            scai::lama::DenseVector<ValueType> dielectricPermittivityEM; //!< Vector storing dielectricPermittivityEM.
            
            scai::lama::DenseVector<ValueType> tauConductivityEM; //!< Vector storing tauConductivityEM for visco-emem modelling.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityEM; //!< Vector storing TauDielectricPermittivityEM for visco-emem modelling.
            scai::IndexType numRelaxationMechanisms; //!< Number of relaxation mechanisms
            ValueType relaxationFrequency;           //!< Relaxation Frequency
            
            scai::lama::DenseVector<ValueType> porosity; //!< Vector storing porosity.
            scai::lama::DenseVector<ValueType> saturation; //!< Vector storing saturation
            
            void initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, ValueType value);
            void initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, std::string filename, scai::IndexType fileFormat);

            void writeParameterisation(scai::lama::Vector<ValueType> const &vector, std::string filename, scai::IndexType fileFormat) const;

            bool normalizeGradient;

          private:
            void allocateParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM);

            void readParameterisation(scai::lama::Vector<ValueType> &vector, std::string filename, scai::IndexType fileFormat);            
        };
    }
}
