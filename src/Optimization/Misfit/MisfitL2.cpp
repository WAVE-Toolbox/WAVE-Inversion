#include "MisfitL2.hpp"

template <typename ValueType>
scai::lama::Scalar KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::Receivers<ValueType> const &receivers1, KITGPI::Acquisition::Receivers<ValueType> const &receivers2)
{    
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandler1;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandler2;
    
    scai::lama::Scalar misfit;

    seismoHandler1 = receivers1.getSeismogramHandler();
    seismoHandler2 = receivers2.getSeismogramHandler();
    misfit = this->calc(seismoHandler1,seismoHandler2);
    
    return misfit;   
}

template <typename ValueType>
scai::lama::Scalar KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::SeismogramHandler<ValueType> const &seismoHandler1, KITGPI::Acquisition::SeismogramHandler<ValueType> const &seismoHandler2)
{
    KITGPI::Acquisition::Seismogram<ValueType> seismogram1;
    KITGPI::Acquisition::Seismogram<ValueType> seismogram2;
    
    scai::lama::Scalar misfit;
    scai::lama::Scalar misfitSum;
    misfitSum = 0;
    
    /* Note that the misfit of different components (p,vx,vy,vz) is summed up. If p and v? is used at the same time, this could cause problems because they could have different scales.
       For different velocity components it should be ok. */
    
    for (int i=0; i<4; i++) {
        seismogram1 = seismoHandler1.getSeismogram(static_cast<KITGPI::Acquisition::SeismogramType>(i));
        seismogram2 = seismoHandler2.getSeismogram(static_cast<KITGPI::Acquisition::SeismogramType>(i));
        misfit = this->calc(seismogram1,seismogram2);
        misfitSum += misfit;
    }
    
    return misfitSum;  
}

template <typename ValueType>
scai::lama::Scalar KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::Seismogram<ValueType> &seismogram1, KITGPI::Acquisition::Seismogram<ValueType> &seismogram2)
{    
    KITGPI::Acquisition::Seismogram<ValueType> seismogramDiff;
    
    seismogramDiff = seismogram1 - seismogram2;
    return 0.5*seismogramDiff.getData().l2Norm();  
}


template class KITGPI::Misfit::MisfitL2<double>;
template class KITGPI::Misfit::MisfitL2<float>;
