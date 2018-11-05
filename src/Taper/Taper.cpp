#include "Taper.hpp"

using namespace scai;

/*! \brief Wrapper to support SeismogramHandler
 \param seismograms SeismogramHandler object
 */
template <typename ValueType>
void KITGPI::Taper::Taper<ValueType>::apply(KITGPI::Acquisition::SeismogramHandler<ValueType> &seismograms) const {
    for (scai::IndexType iComponent = 0; iComponent < 4; iComponent++) {
        if (seismograms.getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            Acquisition::Seismogram<ValueType> &thisSeismogram = seismograms.getSeismogram(Acquisition::SeismogramType(iComponent));
            apply(thisSeismogram);
        }
    }
}

template <typename ValueType>
void KITGPI::Taper::Taper<ValueType>::init(dmemo::DistributionPtr , hmemo::ContextPtr , bool ) {
    COMMON_THROWEXCEPTION("unsupported taper init")
}

template <typename ValueType>
void KITGPI::Taper::Taper<ValueType>::init(dmemo::DistributionPtr ,dmemo::DistributionPtr , hmemo::ContextPtr ) {
    COMMON_THROWEXCEPTION("unsupported taper init")
}

template <typename ValueType>
bool KITGPI::Taper::Taper<ValueType>::getDirection() const {
    COMMON_THROWEXCEPTION("this taper has no direction")
}

template <typename ValueType>
void KITGPI::Taper::Taper<ValueType>::calcCosineTaper(IndexType , IndexType , bool ) {
    COMMON_THROWEXCEPTION("cosine unsupported for this taper type")
}

template <typename ValueType>
void KITGPI::Taper::Taper<ValueType>::calcCosineTaper(IndexType , IndexType , IndexType , IndexType , bool ) {
    COMMON_THROWEXCEPTION("cosine unsupported for this taper type")
}

template <typename ValueType>
void KITGPI::Taper::Taper<ValueType>::calcCosineTaperUp(lama::DenseVector<ValueType> &, IndexType , IndexType ) {
    COMMON_THROWEXCEPTION("cosine unsupported for this taper type")
}

template <typename ValueType>
void KITGPI::Taper::Taper<ValueType>::calcCosineTaperDown(lama::DenseVector<ValueType> &, IndexType , IndexType ) {
    COMMON_THROWEXCEPTION("cosine unsupported for this taper type")
}

template <typename ValueType>
void KITGPI::Taper::Taper<ValueType>::readTaper(std::string , IndexType ) {
    COMMON_THROWEXCEPTION("reading of this taper type not supported yet")
}

template <typename ValueType>
void KITGPI::Taper::Taper<ValueType>::apply(KITGPI::Gradient::Gradient<ValueType> & ) const{
    COMMON_THROWEXCEPTION("this taper can't be applied to a gradient")
}

template class KITGPI::Taper::Taper<double>;
template class KITGPI::Taper::Taper<float>;