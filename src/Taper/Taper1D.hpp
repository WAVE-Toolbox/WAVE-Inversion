#pragma once
#include "Taper.hpp"

namespace KITGPI
{

    namespace Taper
    {

        //! \brief 1-D Taper
        template <typename ValueType>
        class Taper1D : public Taper<ValueType>
        {

          public:
            //! Default constructor
            Taper1D() : direction(0){};

            //! Default destructor
            ~Taper1D(){};

            void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, bool dir) override;

            void calcCosineTaper(scai::IndexType iStart, scai::IndexType iEnd, bool reverse) override;
            void calcCosineTaper(scai::IndexType iStart1, scai::IndexType iEnd1, scai::IndexType iStart2, scai::IndexType iEnd2, bool reverse) override;

            bool getDirection() const override;

            void apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const;
            void apply(KITGPI::Gradient::Gradient<ValueType> &grad) const override;
            void apply(scai::lama::DenseMatrix<ValueType> &mat) const;

            void read(std::string filename, scai::IndexType partitionedIn) override;

          private:
            void calcCosineTaperUp(scai::lama::DenseVector<ValueType> &result, scai::IndexType iStart, scai::IndexType iEnd) override;
            void calcCosineTaperDown(scai::lama::DenseVector<ValueType> &result, scai::IndexType iStart, scai::IndexType iEnd) override;

            bool direction; // 1D taper direction (0=vertical, 1=horizontal)

            scai::lama::DenseVector<ValueType> data;
        };
    }
}