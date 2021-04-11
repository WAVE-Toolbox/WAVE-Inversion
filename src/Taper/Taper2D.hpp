#pragma once
#include <cmath>
#include <algorithm>
#include <vector>     
#include <scai/lama.hpp>
#include <Acquisition/SeismogramHandler.hpp>
#include <Wavefields/WavefieldsFactory.hpp>
#include <Configuration/Configuration.hpp>

namespace KITGPI
{

    namespace Taper
    {

        //! \brief 2-D Taper
        template <typename ValueType>
        class Taper2D
        {

          public:
            //! Default constructor
            Taper2D(){};

            //! Default destructor
            ~Taper2D(){};

            void init(scai::dmemo::DistributionPtr rowDist, scai::dmemo::DistributionPtr colDist, scai::hmemo::ContextPtr ctx);
            void init(KITGPI::Acquisition::SeismogramHandler<ValueType> const seismograms);
            void initWavefieldTransform(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr averageDist, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, bool isSeismic);

            void apply(KITGPI::Acquisition::SeismogramHandler<ValueType> &seismograms) const;
            void apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const;
            void apply(scai::lama::DenseMatrix<ValueType> &mat) const;
            
            KITGPI::Wavefields::Wavefields<ValueType> &applyWavefieldAverage(typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr &wavefieldPtr);
            KITGPI::Wavefields::Wavefields<ValueType> &applyWavefieldRecover(typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr &wavefieldPtr);

            void read(std::string filename);
            void calcInversionAverageMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesInversion);

          private:
            scai::lama::DenseMatrix<ValueType> data;
            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Define sparse format as CSRSparseMatrix
            SparseFormat wavefieldAverageMatrix;
            SparseFormat wavefieldRecoverMatrix;
            typedef typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
            wavefieldPtr wavefieldAverage;
            wavefieldPtr wavefieldRecover;
        };
    }
}
