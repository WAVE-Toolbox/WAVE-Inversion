

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

            /*! \brief Abstract initialisation function
             * Standard initialisation function
             \param ctx Context
             \param dist Distribution
             */
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            virtual void resetGradient() = 0;

            virtual void estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;

            /*! \brief Abstract write function
             *
             * Standard write function
             *
             \param filename filename to write parameters (endings will be added by derived classes)
             \param fileFormat format of output file
             */
            virtual void write(std::string filename, scai::IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const = 0;

            virtual std::string getEquationType() const = 0;

            virtual scai::lama::Vector<ValueType> const &getDensity();
            virtual scai::lama::Vector<ValueType> const &getDensity() const;
            virtual scai::lama::Vector<ValueType> const &getVelocityP();
            virtual scai::lama::Vector<ValueType> const &getVelocityP() const;
            virtual scai::lama::Vector<ValueType> const &getVelocityS();
            virtual scai::lama::Vector<ValueType> const &getVelocityS() const;

            virtual scai::lama::Vector<ValueType> const &getTauP();
            virtual scai::lama::Vector<ValueType> const &getTauP() const;
            virtual scai::lama::Vector<ValueType> const &getTauS();
            virtual scai::lama::Vector<ValueType> const &getTauS() const;
            bool getNormalizeGradient() const;

            virtual scai::IndexType getNumRelaxationMechanisms() const;
            virtual ValueType getRelaxationFrequency() const;

            virtual void setDensity(scai::lama::Vector<ValueType> const &setDensity);
            virtual void setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP);
            virtual void setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS);

            virtual void setTauP(scai::lama::Vector<ValueType> const &setTauP);
            virtual void setTauS(scai::lama::Vector<ValueType> const &setTauS);

            virtual void setNumRelaxationMechanisms(scai::IndexType const setNumRelaxationMechanisms);
            virtual void setRelaxationFrequency(ValueType const setRelaxationFrequency);

            void setNormalizeGradient(bool const &normGrad);

            virtual void scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow) = 0;
            virtual void normalize() = 0;

            virtual void minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs) = 0;
            virtual void plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs) = 0;
            virtual void assign(KITGPI::Gradient::Gradient<ValueType> const &rhs) = 0;
            virtual void timesAssign(ValueType const &rhs) = 0;
            virtual void timesAssign(scai::lama::Vector<ValueType> const &rhs) = 0;

            virtual void sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot) = 0;
            
            virtual void setGradientSubset(KITGPI::Gradient::Gradient<ValueType> &gradientSmall, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoord, scai::IndexType cutCoordInd, scai::IndexType smoothRange, scai::IndexType NX, scai::IndexType NY, scai::IndexType NXBig, scai::IndexType NYBig, scai::IndexType boundaryWidth)=0;
            
            virtual void smoothGradient(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType NX, scai::IndexType NY) = 0;

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

          protected:
            void resetParameter(scai::lama::DenseVector<ValueType> &vector) { vector.setScalar(0); }

            std::string equationType;

            scai::lama::DenseVector<ValueType> density; //!< Vector storing Density.

            scai::lama::DenseVector<ValueType> velocityP; //!< Vector storing P-wave velocity.
            scai::lama::DenseVector<ValueType> velocityS; //!< Vector storing S-wave velocity.

            scai::lama::DenseVector<ValueType> tauP; //!< Vector storing tauP for visco-elastic modelling.
            scai::lama::DenseVector<ValueType> tauS; //!< Vector storing tauS for visco-elastic modelling.

            scai::IndexType numRelaxationMechanisms; //!< Number of relaxation mechanisms
            ValueType relaxationFrequency;           //!< Relaxation Frequency

            void initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value);
            void initParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            void writeParameterisation(scai::lama::Vector<ValueType> const &vector, std::string filename, scai::IndexType fileFormat) const;

            bool normalizeGradient;
            
            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Declare Sparse-Matrix
            SparseFormat getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D &cutCoord);
            
            scai::lama::SparseVector<ValueType> getEraseVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D &cutCoord, scai::IndexType NX, scai::IndexType NYBig, scai::IndexType boundaryWidth);
            
            scai::lama::DenseVector<ValueType> smoothParameter(Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, scai::lama::DenseVector<ValueType> parameter, Acquisition::coordinate3D &cutCoord, scai::IndexType smoothRange, scai::IndexType NX, scai::IndexType NXBig, scai::IndexType NYBig);

          private:
            void allocateParameterisation(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void readParameterisation(scai::lama::Vector<ValueType> &vector, std::string filename, scai::IndexType fileFormat);
        };
    }
}
