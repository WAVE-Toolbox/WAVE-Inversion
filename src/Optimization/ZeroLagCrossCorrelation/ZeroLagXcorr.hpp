

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

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

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            virtual void write(std::string type, IndexType t) = 0;

          protected:
            bool invertForVp = true;
            bool invertForVs = false;
            bool invertForDensity = false;

            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            void writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string type, IndexType t);

            scai::lama::DenseVector<ValueType> VSum; //!< Wavefield for velocity in x
            scai::lama::DenseVector<ValueType> P;    //!< Wavefield
        };
    }
}
