

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

#include "../Common/HostPrint.hpp"
#include <Configuration/Configuration.hpp>
#include <PartitionedInOut/PartitionedInOut.hpp>

namespace KITGPI
{

    //! \brief Parameterisation namespace
    namespace Parameterisation
    {

        //! \brief Abstract class for a single Parameterisation (Subsurface properties)
        /*!
         * As this class is an abstract class, all constructors are protected.
         */
        template <typename ValueType>
        class Parameterisation
        {
          public:
            //! Default constructor.
            Parameterisation() :  numRelaxationMechanisms(0){};

            //! Default destructor.
            ~Parameterisation(){};

            //! \brief Parameterisation pointer
            typedef std::shared_ptr<Parameterisation<ValueType>> ParameterisationPtr;

            /*! \brief Abstract initialization function
             * Standard initialisation function
             \param ctx Context
             \param dist Distribution
             \param filename filename to read parameters (endings will be added by derived classes)
             \param partitionedIn Partitioned input
             */
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn) = 0;

            /*! \brief Abstract initialisation function
             * initialises parameter vectors with zero
             \param ctx Context
             \param dist Distribution
             */
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;
	    
            /*! \brief Abstract initialisation function
             * Standard initialisation function
             \param config Configuration from configuration file
             \param ctx Context
             \param dist Distribution
             */
            virtual void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            /*! \brief Abstract write function
             *
             * Standard write function
             *
             \param filename filename to write parameters (endings will be added by derived classes)
             \param partitionedOut Partitioned output
             */
            virtual void write(std::string filename, IndexType partitionedOut) const = 0;

            virtual scai::lama::Vector const &getDensity();
            virtual scai::lama::Vector const &getDensity() const;
            virtual scai::lama::Vector const &getVelocityP();
            virtual scai::lama::Vector const &getVelocityP() const;
            virtual scai::lama::Vector const &getVelocityS();
            virtual scai::lama::Vector const &getVelocityS() const;

            virtual scai::lama::Vector const &getTauP();
            virtual scai::lama::Vector const &getTauP() const;
            virtual scai::lama::Vector const &getTauS();
            virtual scai::lama::Vector const &getTauS() const;

            virtual IndexType getNumRelaxationMechanisms() const;
            virtual ValueType getRelaxationFrequency() const;

            virtual void setDensity(scai::lama::Vector const &setDensity);
            virtual void setVelocityP(scai::lama::Vector const &setVelocityP);
            virtual void setVelocityS(scai::lama::Vector const &setVelocityS);

            virtual void setTauP(scai::lama::Vector const &setTauP);
            virtual void setTauS(scai::lama::Vector const &setTauS);

            virtual void setNumRelaxationMechanisms(IndexType const setNumRelaxationMechanisms);
            virtual void setRelaxationFrequency(ValueType const setRelaxationFrequency);

            virtual void minusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs)=0;
	    virtual void plusAssign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs)=0;
	    virtual void assign(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs)=0;
	    
	    
	    KITGPI::Parameterisation::Parameterisation<ValueType> &operator=(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
	    KITGPI::Parameterisation::Parameterisation<ValueType> &operator-=(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
	    KITGPI::Parameterisation::Parameterisation<ValueType> &operator+=(KITGPI::Parameterisation::Parameterisation<ValueType> const &rhs);
	    
          protected:

            IndexType PartitionedIn;  //!< ==1 If Modulus is read from partitioned fileblock; ==0 if modulus is in single files
            IndexType PartitionedOut; //!< ==1 If Modulus is written to partitioned fileblock; ==0 if modulus is written to single files

            scai::lama::DenseVector<ValueType> density;        //!< Vector storing Density.

            scai::lama::DenseVector<ValueType> velocityP; //!< Vector storing P-wave velocity.
            scai::lama::DenseVector<ValueType> velocityS; //!< Vector storing S-wave velocity.

            scai::lama::DenseVector<ValueType> tauP; //!< Vector storing tauP for visco-elastic modelling.
            scai::lama::DenseVector<ValueType> tauS; //!< Vector storing tauS for visco-elastic modelling.


            IndexType numRelaxationMechanisms; //!< Number of relaxation mechanisms
            ValueType relaxationFrequency;     //!< Relaxation Frequency

            void initParameterisation(scai::lama::Vector &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar value);
            void initParameterisation(scai::lama::Vector &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn);

            void writeParameterisation(scai::lama::Vector const &vector, std::string filename, IndexType partitionedOut) const;


            IndexType getPartitionedIn();
            IndexType getPartitionedOut();

         
   
          private:
            void allocateParameterisation(scai::lama::Vector &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void readParameterisation(scai::lama::Vector &vector, std::string filename, scai::dmemo::DistributionPtr dist, IndexType partitionedIn);


        };
    }
}
