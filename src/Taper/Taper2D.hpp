#pragma once
#include <cmath>
#include <algorithm>
#include <vector>     
#include <scai/lama.hpp>
#include <Acquisition/SeismogramHandler.hpp>
#include <Configuration/Configuration.hpp>
#include <Modelparameter/Modelparameter.hpp>
#include <scai/dmemo/SingleDistribution.hpp>

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
            void initTransformMatrix(scai::dmemo::DistributionPtr dist1, scai::dmemo::DistributionPtr dist2, scai::hmemo::ContextPtr ctx);

            void apply(KITGPI::Acquisition::SeismogramHandler<ValueType> &seismograms) const;
            void apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const;
            void apply(scai::lama::DenseMatrix<ValueType> &mat) const;
            
            void applyGradientTransform1to2(scai::lama::Vector<ValueType> const &gradientParameter1, lama::Vector<ValueType> &gradientParameter2);
            void applyGradientTransform2to1(scai::lama::Vector<ValueType> &gradientParameter1, scai::lama::Vector<ValueType> const &gradientParameter2);
            void applyModelTransform1to2(scai::lama::Vector<ValueType> const &modelParameter1, scai::lama::Vector<ValueType> &modelParameter2);
            void applyModelTransform2to1(scai::lama::Vector<ValueType> &modelParameter1, scai::lama::Vector<ValueType> const &modelParameter2);

            void read(std::string filename);
            void calcTransformMatrix1to2(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates1, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates2);
            void calcTransformMatrix2to1(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates1, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates2);

            void exchangePetrophysics(KITGPI::Modelparameter::Modelparameter<ValueType> const &model2, KITGPI::Modelparameter::Modelparameter<ValueType> &model1, KITGPI::Configuration::Configuration config, IndexType equationInd);
            void exchangeModelparameters(KITGPI::Modelparameter::Modelparameter<ValueType> const &model2, KITGPI::Configuration::Configuration config2, KITGPI::Modelparameter::Modelparameter<ValueType> &model1, KITGPI::Configuration::Configuration config1, IndexType equationInd);
            
            void initAverageMatrix(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr distInversion, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx);
            void calcAverageMatrix(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates, KITGPI::Acquisition::Coordinates<ValueType> modelCoordinatesInversion);
            scai::lama::Matrix<ValueType> const &getAverageMatrix();
            scai::lama::Matrix<ValueType> const &getRecoverMatrix();
            
            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Define sparse format as CSRSparseMatrix
            SparseFormat averageMatrix;
            SparseFormat recoverMatrix;
            
          private:
            scai::lama::DenseMatrix<ValueType> data;
            SparseFormat transformMatrix1to2;
            SparseFormat transformMatrix2to1;            
        };
    }
}
