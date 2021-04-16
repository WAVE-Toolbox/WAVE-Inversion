#pragma once
#include <cmath>
#include <algorithm>
#include <vector>     
#include <scai/lama.hpp>
#include <Acquisition/SeismogramHandler.hpp>
#include <AcquisitionEM/SeismogramHandler.hpp>
#include <Wavefields/WavefieldsFactory.hpp>
#include <WavefieldsEM/WavefieldsFactory.hpp>
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
            void init(KITGPI::Acquisition::SeismogramHandlerEM<ValueType> const seismograms);
            void initModelTransform(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distEM, scai::hmemo::ContextPtr ctx);
            void initWavefieldTransform(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr distInversion, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, bool isSeismic);

            void apply(KITGPI::Acquisition::SeismogramHandler<ValueType> &seismograms) const;
            void apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const;
            void apply(KITGPI::Acquisition::SeismogramHandlerEM<ValueType> &seismograms) const;
            void apply(KITGPI::Acquisition::SeismogramEM<ValueType> &seismogram) const;
            void apply(scai::lama::DenseMatrix<ValueType> &mat) const;
            
            scai::lama::Vector<ValueType> const &applyGradientTransformToEM(scai::lama::Vector<ValueType> const &gradientParameter);
            scai::lama::Vector<ValueType> const &applyGradientTransformToSeismic(scai::lama::Vector<ValueType> const &gradientParameterEM);
            scai::lama::Vector<ValueType> const &applyModelTransformToEM(scai::lama::Vector<ValueType> const &modelParameter, scai::lama::Vector<ValueType> const &modelParameterEM);
            scai::lama::Vector<ValueType> const &applyModelTransformToSeismic(scai::lama::Vector<ValueType> const &modelParameter, scai::lama::Vector<ValueType> const &modelParameterEM);
            KITGPI::Wavefields::Wavefields<ValueType> &applyWavefieldAverage(typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr &wavefieldPtr);
            KITGPI::Wavefields::WavefieldsEM<ValueType> &applyWavefieldAverage(typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr &wavefieldPtrEM);
            KITGPI::Wavefields::Wavefields<ValueType> &applyWavefieldRecover(typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr &wavefieldPtr);
            KITGPI::Wavefields::WavefieldsEM<ValueType> &applyWavefieldRecover(typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr &wavefieldPtrEM);

            void read(std::string filename);
            void calcInversionAverageMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesInversion);
            void calcSeismictoEMMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesEM, KITGPI::Configuration::Configuration configEM);
            void calcEMtoSeismicMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesEM, KITGPI::Configuration::Configuration configEM);

            void exchangePetrophysics(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelEM, KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Configuration::Configuration config);
            void exchangePetrophysics(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelEM, KITGPI::Configuration::Configuration configEM);
            
          private:
            scai::lama::DenseMatrix<ValueType> data;
            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Define sparse format as CSRSparseMatrix
            SparseFormat modelTransformMatrixToEM;
            SparseFormat modelTransformMatrixToSeismic;
            SparseFormat wavefieldAverageMatrix;
            SparseFormat wavefieldRecoverMatrix;
            scai::lama::DenseVector<ValueType> modelParameterTransform;
            scai::lama::DenseVector<ValueType> modelParameterTransformEM;
            typedef typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
            wavefieldPtr wavefieldAverage;
            wavefieldPtr wavefieldRecover;
            typedef typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr wavefieldPtrEM; 
            wavefieldPtrEM wavefieldAverageEM;
            wavefieldPtrEM wavefieldRecoverEM;
        };
    }
}
