#pragma once
#include <cmath>
#include <algorithm>
#include <vector>     
#include <scai/lama.hpp>
#include <Acquisition/SeismogramHandler.hpp>
#include <Configuration/Configuration.hpp>
#include <Modelparameter/Modelparameter.hpp>
#include <Wavefields/WavefieldsFactory.hpp>

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
            void initModelTransform(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distEM, scai::hmemo::ContextPtr ctx);

            void apply(KITGPI::Acquisition::SeismogramHandler<ValueType> &seismograms) const;
            void apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const;
            void apply(scai::lama::DenseMatrix<ValueType> &mat) const;
            
            scai::lama::Vector<ValueType> const &applyGradientTransformToEM(scai::lama::Vector<ValueType> const &gradientParameter);
            scai::lama::Vector<ValueType> const &applyGradientTransformToSeismic(scai::lama::Vector<ValueType> const &gradientParameterEM);
            scai::lama::Vector<ValueType> const &applyModelTransformToEM(scai::lama::Vector<ValueType> const &modelParameter, scai::lama::Vector<ValueType> const &modelParameterEM);
            scai::lama::Vector<ValueType> const &applyModelTransformToSeismic(scai::lama::Vector<ValueType> const &modelParameter, scai::lama::Vector<ValueType> const &modelParameterEM);

            void read(std::string filename);
            void calcSeismictoEMMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesEM, KITGPI::Configuration::Configuration configEM);
            void calcEMtoSeismicMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesEM, KITGPI::Configuration::Configuration configEM);

            void exchangePetrophysics(KITGPI::Modelparameter::Modelparameter<ValueType> &model, KITGPI::Modelparameter::Modelparameter<ValueType> &modelEM, KITGPI::Configuration::Configuration config);
            void exchangeModelparameters(KITGPI::Modelparameter::Modelparameter<ValueType> &model1, KITGPI::Modelparameter::Modelparameter<ValueType> &model2, KITGPI::Configuration::Configuration config1, KITGPI::Configuration::Configuration config2);
            
            void initWavefieldAverageMatrix(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr distInversion, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx);
            void calcWavefieldAverageMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesInversion);
            KITGPI::Wavefields::Wavefields<ValueType> &applyWavefieldAverage(typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr &wavefieldPtr);
            KITGPI::Wavefields::Wavefields<ValueType> &applyWavefieldRecover(typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr &wavefieldPtr);
            
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
            wavefieldPtr wavefieldsInversion;
            wavefieldPtr wavefields;
        };
    }
}
