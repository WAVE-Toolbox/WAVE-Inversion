#include "Taper2D.hpp"

using namespace scai;

/*! \brief Initialize taper
 \param rowDist Row distribution
 \param colDist Column distribution
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::init(dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist, hmemo::ContextPtr ctx) {
    
    data.allocate(rowDist, colDist);
    data.setContextPtr(ctx);
}

/*! \brief Apply taper to a single seismogram
 \param seismogram Seismogram
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const {
    lama::DenseMatrix<ValueType> &seismogramData = seismogram.getData();
    seismogramData.binaryOp(seismogramData,common::BinaryOp::MULT,data);
}

/*! \brief Apply taper to a DenseMatrix
 \param seismogram Seismogram
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::apply(lama::DenseMatrix<ValueType> &mat) const {
    mat.binaryOp(mat,common::BinaryOp::MULT,data);
}

/*! \brief Read a taper from file
 * \param filename taper filename
 */
template <typename ValueType>
void KITGPI::Taper::Taper2D<ValueType>::readTaper(std::string filename, IndexType /*partitionedIn*/)
{
    scai::dmemo::DistributionPtr distTraces(data.getRowDistributionPtr());
    scai::dmemo::DistributionPtr distSamples(data.getColDistributionPtr());
    
    data.readFromFile(filename);
    
    data.redistribute(distTraces, distSamples);
}

template class KITGPI::Taper::Taper2D<double>;
template class KITGPI::Taper::Taper2D<float>;