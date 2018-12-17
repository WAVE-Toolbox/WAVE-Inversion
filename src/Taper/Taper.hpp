#pragma once
#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include <Configuration/Configuration.hpp>
#include <Acquisition/SeismogramHandler.hpp>
#include <PartitionedInOut/PartitionedInOut.hpp>
#include "../Gradient/Gradient.hpp"

namespace KITGPI
{

    //! \brief Taper namespace
    namespace Taper
    {

        //! \brief Abstract class for Taper
        template <typename ValueType>
        class Taper
        {
            public:
                //! \brief Declare Taper pointer
                typedef std::shared_ptr<Taper<ValueType>> TaperPtr;
            
                //! Default constructor
                Taper(){};
                
                //! Default destructor
                ~Taper(){};
                
                virtual void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, bool dir);
                virtual void init(scai::dmemo::DistributionPtr rowDist, scai::dmemo::DistributionPtr colDist, scai::hmemo::ContextPtr ctx);
                
                virtual void calcCosineTaper(scai::IndexType iStart, scai::IndexType iEnd, bool reverse);
                virtual void calcCosineTaper(scai::IndexType iStart1, scai::IndexType iEnd1, scai::IndexType iStart2, scai::IndexType iEnd2, bool reverse);
                
                virtual bool getDirection() const;
                
                void apply(KITGPI::Acquisition::SeismogramHandler<ValueType> &seismograms) const;
                virtual void apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const = 0;
                virtual void apply(KITGPI::Gradient::Gradient<ValueType> &grad) const;
                virtual void apply(scai::lama::DenseMatrix<ValueType> &mat) const = 0;
                
                virtual void readTaper(std::string filename, scai::IndexType partitionedIn) = 0;
                
            private:
                
                virtual void calcCosineTaperUp(scai::lama::DenseVector<ValueType> &result, scai::IndexType iStart, scai::IndexType iEnd);
                virtual void calcCosineTaperDown(scai::lama::DenseVector<ValueType> &result, scai::IndexType iStart, scai::IndexType iEnd); 
              
        };
    }
}