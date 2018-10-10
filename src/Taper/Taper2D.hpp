#pragma once
#include "Taper.hpp"

namespace KITGPI
{

    namespace Taper
    {

        //! \brief 2-D Taper
        template <typename ValueType>
        class Taper2D : public Taper<ValueType>
        {

        public:
            //! Default constructor
            Taper2D(){};

            //! Default destructor
            ~Taper2D(){};
            
            void init(scai::dmemo::DistributionPtr rowDist, scai::dmemo::DistributionPtr colDist, scai::hmemo::ContextPtr ctx) override;
            
            void apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const;
            
        private:
            
            scai::lama::DenseMatrix<ValueType> data;
            
        };
    }
}