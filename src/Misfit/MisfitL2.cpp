#include "MisfitL2.hpp"
#include <Acquisition/Acquisition.hpp>

/*! \brief Return the L2-norm of seismograms stored in the given receiver objects (note: the misfit is summed over all components, i.e. vx,vy,vz,p)
 *
 *
 \param receiversSyn Receiver object which stores the synthetic data 
 \param receiversObs Receiver object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs)
{    
        
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObs;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerSyn;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerObs;
    
    ValueType misfit;
    ValueType misfitSum;
    misfitSum = 0;

    seismoHandlerSyn = receiversSyn.getSeismogramHandler();
    seismoHandlerObs = receiversObs.getSeismogramHandler();
      
    /* Note that the misfit of different components (p,vx,vy,vz) is summed up. If p and v? is used at the same time, this could cause problems because they could have different scales.
       For different velocity components it should be ok. */
    
    for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<KITGPI::Acquisition::SeismogramType>(i));
        seismogramObs = seismoHandlerObs.getSeismogram(static_cast<KITGPI::Acquisition::SeismogramType>(i));
        misfit = this->calc(seismogramSyn, seismogramObs);
        misfitSum += misfit;
    }
    
    return misfitSum; 
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");
    
    KITGPI::Acquisition::Seismogram<ValueType> seismogramDiff;
    
    seismogramDiff = seismogramSyn - seismogramObs;
    return 0.5*seismogramDiff.getData().l2Norm();  
}

/*! \brief Calculate the adjoint sources
 *
 *
 \param adjointSources Receiver object which stores the adjoint sources (in- and output)
 \param receiversSyn Receiver object which stores the synthetic data 
 \param receiversObs Receiver object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSources(KITGPI::Acquisition::Receivers<ValueType> &adjointSources, KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs)
{
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObs;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramAdj;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerSyn;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerObs;
    
    seismoHandlerSyn = receiversSyn.getSeismogramHandler();
    seismoHandlerObs = receiversObs.getSeismogramHandler();
    
    for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<KITGPI::Acquisition::SeismogramType>(i));
        seismogramObs = seismoHandlerObs.getSeismogram(static_cast<KITGPI::Acquisition::SeismogramType>(i));
        this->calcAdjointSeismogram(seismogramAdj, seismogramSyn, seismogramObs);
        adjointSources.getSeismogramHandler().getSeismogram(seismogramAdj.getTraceType()) = seismogramAdj;
    }    
    
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogram(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{
    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    seismogramAdj = seismogramSyn - seismogramObs; 

}


template class KITGPI::Misfit::MisfitL2<double>;
template class KITGPI::Misfit::MisfitL2<float>;
