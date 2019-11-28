#pragma once
#include <scai/lama.hpp>
#include <Acquisition/Seismogram.hpp>
#include <Acquisition/SeismogramHandler.hpp>
#include "../Gradient/Gradient.hpp"

namespace KITGPI
{

    namespace Taper
    {

        //! \brief 1-D Taper
        template <typename ValueType>
        class Taper1D
        {

          public:
            //! Default constructor
            Taper1D(){};

            //! Default destructor
            ~Taper1D(){};

            void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, bool dir);

            void calcCosineTaper(scai::IndexType iStart, scai::IndexType iEnd, bool reverse);
            void calcCosineTaper(scai::IndexType iStart1, scai::IndexType iEnd1, scai::IndexType iStart2, scai::IndexType iEnd2, bool reverse);

            void apply(KITGPI::Acquisition::SeismogramHandler<ValueType> &seismograms) const;
            void apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const;
            void apply(KITGPI::Gradient::Gradient<ValueType> &grad) const;
            void apply(scai::lama::DenseMatrix<ValueType> &mat) const;

            void read(std::string filename, scai::IndexType fileFormat);
            
            bool getDirection() const;

          private:
            void calcCosineTaperUp(scai::lama::DenseVector<ValueType> &result, scai::IndexType iStart, scai::IndexType iEnd);
            void calcCosineTaperDown(scai::lama::DenseVector<ValueType> &result, scai::IndexType iStart, scai::IndexType iEnd);

            scai::lama::DenseVector<ValueType> data;
            bool direction;
        };
    }
}
