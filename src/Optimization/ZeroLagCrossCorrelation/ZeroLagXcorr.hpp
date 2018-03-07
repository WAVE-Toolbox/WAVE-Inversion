

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <Configuration/Configuration.hpp>
#include <Wavefields/WavefieldsFactory.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

namespace KITGPI
{

    //! \brief ZeroLagXcorr namespace
    namespace ZeroLagXcorr
    {

        /*! \brief Abstract class to calculate the cross correlation between Forward and Adjoint Wavefield.
         *
         * ZeroLagXcorr implements some methods, which are requiered by all derived classes.
         * As this class is an abstract class, all methods are protected.
         */
        template <typename ValueType>
        class ZeroLagXcorr
        {

          public:
            //! Default constructor
            ZeroLagXcorr(){};
            //! Default destructor
            ~ZeroLagXcorr(){};

            //! \brief Declare ZeroLagXcorr pointer
            typedef std::shared_ptr<ZeroLagXcorr<ValueType>> ZeroLagXcorrPtr;

            //! Reset cross correlated wavefields
            virtual void resetXcorr() = 0;

            virtual void update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield) = 0;

            virtual scai::lama::DenseVector<ValueType> const &getVSum() const;
            virtual scai::lama::DenseVector<ValueType> const &getP() const;
            virtual scai::lama::DenseVector<ValueType> const &getShearStress() const;
            virtual scai::lama::DenseVector<ValueType> const &getNormalStressDiff() const;
            virtual scai::lama::DenseVector<ValueType> const &getNormalStressSum() const;

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            virtual void write(std::string type, IndexType t) = 0;

            bool invertForVp = false;
            bool invertForVs = false;
            bool invertForDensity = false;

          protected:
            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            void writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string type, IndexType t);

            scai::lama::DenseVector<ValueType> ShearStress;      //!< (sum of) correlated shear stresses
            scai::lama::DenseVector<ValueType> NormalStressDiff; //!<correlated difference of normal stress components  2D: (sxxF-syyF)*(sxxB-syyB)
            scai::lama::DenseVector<ValueType> NormalStressSum;  //!<correlated sum        of normal stress components  2D: (sxxF+syyF)*(sxxB+syyB)
            scai::lama::DenseVector<ValueType> VSum;             //!< sum of the correlated velocity wavefields sum_i (viF*viB)
            scai::lama::DenseVector<ValueType> P;                //!< correlated pressure Wavefield
        };
    }
}
