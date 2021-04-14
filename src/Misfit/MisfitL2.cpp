#include "MisfitL2.hpp"

/*! \brief Return the L2-norm of seismograms stored in the given receiver objects (note: the misfit is summed over all components, i.e. vx,vy,vz,p)
 *
 *
 \param receiversSyn Receiver object which stores the synthetic data 
 \param receiversObs Receiver object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs, scai::IndexType shotInd)
{            
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObs;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerSyn;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerObs;
    
    ValueType misfit = 0;
    ValueType misfitSum = 0;
    scai::IndexType count = 0;

    seismoHandlerSyn = receiversSyn.getSeismogramHandler();
    seismoHandlerObs = receiversObs.getSeismogramHandler();
              
    /* Note that the misfit of different components (p,vx,vy,vz) is summed up. If p and v? is used at the same time, this could cause problems because they could have different scales.
       For different velocity components it should be ok. */
    
    std::string misfitTypeL2 = this->getMisfitType(shotInd);
    for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
        seismogramObs = seismoHandlerObs.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
        if (misfitTypeL2.compare("l2") == 0) {
            misfit = this->calcL2(seismogramSyn, seismogramObs);
        } else if (misfitTypeL2.compare("l7") == 0) {
            misfit = this->calcL2Normalized(seismogramSyn, seismogramObs);
        } else if (misfitTypeL2.compare("l8") == 0) {
            misfit = this->calcL2Envelope(seismogramSyn, seismogramObs);
        } 
        misfitSum += misfit;
        if (misfit != 0) count++;
    }
    
    return misfitSum/count/this->getMisfitTypeHistory().size(); 
}

/*! \brief Calculate the adjoint sources
 *
 *
 \param adjointSources Receiver object which stores the adjoint sources (in- and output)
 \param receiversSyn Receiver object which stores the synthetic data 
 \param receiversObs Receiver object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSources(KITGPI::Acquisition::Receivers<ValueType> &adjointSources, KITGPI::Acquisition::Receivers<ValueType> const &receiversSyn, KITGPI::Acquisition::Receivers<ValueType> const &receiversObs, scai::IndexType shotInd)
{
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObs;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramAdj;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerSyn;
    KITGPI::Acquisition::SeismogramHandler<ValueType> seismoHandlerObs;
    
    seismoHandlerSyn = receiversSyn.getSeismogramHandler();
    seismoHandlerObs = receiversObs.getSeismogramHandler();
    
    std::string misfitTypeL2 = this->getMisfitType(shotInd);
    for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
        seismogramObs = seismoHandlerObs.getSeismogram(static_cast<Acquisition::SeismogramType>(i));
        if (misfitTypeL2.compare("l2") == 0) {
            this->calcAdjointSeismogramL2(seismogramAdj, seismogramSyn, seismogramObs);
        } else if (misfitTypeL2.compare("l7") == 0) {
            this->calcAdjointSeismogramL2Normalized(seismogramAdj, seismogramSyn, seismogramObs);
        } else if (misfitTypeL2.compare("l8") == 0) {
            this->calcAdjointSeismogramL2Envelope(seismogramAdj, seismogramSyn, seismogramObs);
        } 
        adjointSources.getSeismogramHandler().getSeismogram(seismogramAdj.getTraceType()) = seismogramAdj;
    }       
}

/*! \brief Return the L2-norm of seismograms stored in the given receiver objects (note: the misfit is summed over all components, i.e. vx,vy,vz,p)
 *
 *
 \param receiversSyn Receiver object which stores the synthetic data 
 \param receiversObs Receiver object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calc(KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversSyn, KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversObs, scai::IndexType shotInd)
{    
        
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObs;
    KITGPI::Acquisition::SeismogramHandlerEM<ValueType> seismoHandlerSyn;
    KITGPI::Acquisition::SeismogramHandlerEM<ValueType> seismoHandlerObs;
    
    ValueType misfit = 0;
    ValueType misfitSum = 0;
    scai::IndexType count = 0;

    seismoHandlerSyn = receiversSyn.getSeismogramHandler();
    seismoHandlerObs = receiversObs.getSeismogramHandler();
          
    /* Note that the misfit of different components (p,vx,vy,vz) is summed up. If p and v? is used at the same time, this could cause problems because they could have different scales.
       For different velocity components it should be ok. */
    
    std::string misfitTypeL2 = this->getMisfitType(shotInd);
    for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<Acquisition::SeismogramTypeEM>(i));
        seismogramObs = seismoHandlerObs.getSeismogram(static_cast<Acquisition::SeismogramTypeEM>(i));
        if (misfitTypeL2.compare("l2") == 0) {
            misfit = this->calcL2(seismogramSyn, seismogramObs);
        } else if (misfitTypeL2.compare("l7") == 0) {
            misfit = this->calcL2Normalized(seismogramSyn, seismogramObs);
        } else if (misfitTypeL2.compare("l8") == 0) {
            misfit = this->calcL2Envelope(seismogramSyn, seismogramObs);
        } 
        misfitSum += misfit;
        if (misfit != 0) count++;
    }
    
    return misfitSum/count/this->getMisfitTypeHistory().size(); 
}

/*! \brief Calculate the adjoint sourcesEM
 *
 *
 \param adjointSourcesEM Receiver object which stores the adjoint sourcesEM (in- and output)
 \param receiversSyn Receiver object which stores the synthetic data 
 \param receiversObs Receiver object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSources(KITGPI::Acquisition::ReceiversEM<ValueType> &adjointSourcesEM, KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversSyn, KITGPI::Acquisition::ReceiversEM<ValueType> const &receiversObs, scai::IndexType shotInd)
{
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObs;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramAdj;
    KITGPI::Acquisition::SeismogramHandlerEM<ValueType> seismoHandlerSyn;
    KITGPI::Acquisition::SeismogramHandlerEM<ValueType> seismoHandlerObs;
    
    seismoHandlerSyn = receiversSyn.getSeismogramHandler();
    seismoHandlerObs = receiversObs.getSeismogramHandler();
    
    std::string misfitTypeL2 = this->getMisfitType(shotInd);
    for (int i=0; i<KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        seismogramSyn = seismoHandlerSyn.getSeismogram(static_cast<Acquisition::SeismogramTypeEM>(i));
        seismogramObs = seismoHandlerObs.getSeismogram(static_cast<Acquisition::SeismogramTypeEM>(i));
        if (misfitTypeL2.compare("l2") == 0) {
            this->calcAdjointSeismogramL2(seismogramAdj, seismogramSyn, seismogramObs);
        } else if (misfitTypeL2.compare("l7") == 0) {
            this->calcAdjointSeismogramL2Normalized(seismogramAdj, seismogramSyn, seismogramObs);
        } else if (misfitTypeL2.compare("l8") == 0) {
            this->calcAdjointSeismogramL2Envelope(seismogramAdj, seismogramSyn, seismogramObs);
        } 
        adjointSourcesEM.getSeismogramHandler().getSeismogram(seismogramAdj.getTraceType()) = seismogramAdj;
    }       
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    
    if (seismogramSyn.getData().getNumRows()!=0) {
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows(); 
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm; 
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    seismogramAdj = seismogramSyn - seismogramObs; 
    
    if (seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) {
        seismogramAdj *= -1;        
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2(KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    
    if (seismogramSyn.getData().getNumRows()!=0) {
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows();  
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;    
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2(KITGPI::Acquisition::SeismogramEM<ValueType> &seismogramAdj, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    seismogramAdj = seismogramSyn - seismogramObs; 
    
    if (seismogramSyn.getTraceType() != Acquisition::SeismogramTypeEM::HZ) {
        seismogramAdj *= -1;        
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2Normalized(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    // see Groos L. 2D full waveform inversion of shallow seismic Rayleigh waves[D]. Verlag nicht ermittelbar, 2013 for details
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramSyntemp.normalizeTraceL2();
        seismogramObstemp.normalizeTraceL2();
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows(); 
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;   
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2Normalized(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    // see Groos L. 2D full waveform inversion of shallow seismic Rayleigh waves[D]. Verlag nicht ermittelbar, 2013 for details
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramAdj = seismogramSyn;
        seismogramAdj *= seismogramObs;
        scai::lama::DenseVector<ValueType> tempL2NormSyn = seismogramSyntemp.getTraceL2norm();
        scai::lama::DenseVector<ValueType> tempL2NormObs = seismogramObstemp.getTraceL2norm();
        scai::lama::DenseVector<ValueType> tempSumSynObs = seismogramAdj.getTraceSum();
        tempL2NormSyn = 1 / tempL2NormSyn;
        tempL2NormObs = 1 / tempL2NormObs;     
        tempL2NormObs *= tempL2NormSyn;
        tempL2NormObs *= tempL2NormSyn;
        tempL2NormObs *= tempSumSynObs;
        seismogramSyntemp.normalizeTraceL2();
        seismogramObstemp.normalizeTraceL2();
        seismogramSyntemp.getData().scaleRows(tempL2NormObs);
        seismogramObstemp.getData().scaleRows(tempL2NormSyn);      
    }
    seismogramAdj = seismogramSyntemp - seismogramObstemp;
    
    if (seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) {
        seismogramAdj *= -1;        
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2Normalized(KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramSyntemp.normalizeTraceL2();
        seismogramObstemp.normalizeTraceL2();
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows(); 
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;   
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2Normalized(KITGPI::Acquisition::SeismogramEM<ValueType> &seismogramAdj, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramAdj = seismogramSyn;
        seismogramAdj *= seismogramObs;
        scai::lama::DenseVector<ValueType> tempL2NormSyn = seismogramSyntemp.getTraceL2norm();
        scai::lama::DenseVector<ValueType> tempL2NormObs = seismogramObstemp.getTraceL2norm();
        scai::lama::DenseVector<ValueType> tempSumSynObs = seismogramAdj.getTraceSum();
        tempL2NormSyn = 1 / tempL2NormSyn;
        tempL2NormObs = 1 / tempL2NormObs;     
        tempL2NormObs *= tempL2NormSyn;
        tempL2NormObs *= tempL2NormSyn;
        tempL2NormObs *= tempSumSynObs;
        seismogramSyntemp.normalizeTraceL2();
        seismogramObstemp.normalizeTraceL2();
        seismogramSyntemp.getData().scaleRows(tempL2NormObs);
        seismogramObstemp.getData().scaleRows(tempL2NormSyn);      
    }
    seismogramAdj = seismogramSyntemp - seismogramObstemp;       
    
    if (seismogramSyn.getTraceType() != Acquisition::SeismogramTypeEM::HZ) {
        seismogramAdj *= -1;        
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2Envelope(KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    // see Yuan Y O, Simons F J, Bozdağ E. Multiscale adjoint waveform tomography for surface and body waves[J]. Geophysics, 2015, 80(5): R281-R302 for details
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        Common::envelope(seismogramSyntemp.getData());
        Common::envelope(seismogramObstemp.getData());
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows(); 
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;   
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2Envelope(KITGPI::Acquisition::Seismogram<ValueType> &seismogramAdj, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramSyn, KITGPI::Acquisition::Seismogram<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    // see Yuan Y O, Simons F J, Bozdağ E. Multiscale adjoint waveform tomography for surface and body waves[J]. Geophysics, 2015, 80(5): R281-R302 for details
    KITGPI::Acquisition::Seismogram<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramAdj = seismogramSyn;
        scai::lama::DenseMatrix<ValueType> tempDataSyn = seismogramSyntemp.getData();
        scai::lama::DenseMatrix<ValueType> tempDataObs = seismogramObstemp.getData();
        scai::lama::DenseVector<ValueType> dataTrace; 
        ValueType waterLevel = 1.0;
        Common::envelope(tempDataSyn);
        Common::envelope(tempDataObs);
        for (int i=0; i<tempDataObs.getNumRows(); i++) {
            tempDataObs.getRow(dataTrace, i);  
            if (i==0) {
                waterLevel = dataTrace.maxNorm();                
            } else {
                if (waterLevel > dataTrace.maxNorm()) waterLevel = dataTrace.maxNorm(); 
            }
        }
        waterLevel *= 1e-1;
        tempDataObs.binaryOp(tempDataSyn, scai::common::BinaryOp::SUB, tempDataObs);
        for (int i=0; i<tempDataSyn.getNumRows(); i++) {
            tempDataSyn.getRow(dataTrace, i);   
            dataTrace += waterLevel;
            tempDataSyn.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
        }
        tempDataObs.binaryOp(tempDataObs, scai::common::BinaryOp::DIVIDE, tempDataSyn);
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, tempDataObs);
        
        tempDataSyn = seismogramAdj.getData();
        Common::hilbert(tempDataSyn);
        tempDataSyn.binaryOp(tempDataSyn, scai::common::BinaryOp::MULT, tempDataObs);
        Common::hilbert(tempDataSyn);
        seismogramObstemp.getData() = tempDataSyn;  
    }
    seismogramAdj = seismogramSyntemp - seismogramObstemp;
    
    if (seismogramSyn.getTraceType() == Acquisition::SeismogramType::P) {
        seismogramAdj *= -1;        
    }
}

/*! \brief Return the L2-norm of seismograms 
 *
 *
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
ValueType KITGPI::Misfit::MisfitL2<ValueType>::calcL2Envelope(KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{    
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObstemp = seismogramObs;
    ValueType tempL2Norm = 0;
    
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        Common::envelope(seismogramSyntemp.getData());
        Common::envelope(seismogramObstemp.getData());
        seismogramSyntemp -= seismogramObstemp;
        
        tempL2Norm = 0.5*seismogramSyntemp.getData().l2Norm() / seismogramSyntemp.getData().getNumColumns() / seismogramSyntemp.getData().getNumRows(); 
        if (isfinite(tempL2Norm)==false) tempL2Norm=0;
    }
    
    return tempL2Norm;   
}

/*! \brief Calculate the adjoint seismograms
 *
 *
 \param seismogramAdj Seismogram object which stores the adjoint data  (in- and output)
 \param seismogramSyn Seismogram object which stores the synthetic data 
 \param seismogramObs Seismogram object which stores the observed data 
 */
template <typename ValueType>
void KITGPI::Misfit::MisfitL2<ValueType>::calcAdjointSeismogramL2Envelope(KITGPI::Acquisition::SeismogramEM<ValueType> &seismogramAdj, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramSyn, KITGPI::Acquisition::SeismogramEM<ValueType> const &seismogramObs)
{
    SCAI_ASSERT_ERROR(seismogramSyn.getTraceType() == seismogramObs.getTraceType(), "Seismogram types differ!");

    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramSyntemp = seismogramSyn;
    KITGPI::Acquisition::SeismogramEM<ValueType> seismogramObstemp = seismogramObs;
    if (seismogramSyntemp.getData().getNumRows()!=0) {
        seismogramAdj = seismogramSyn;
        scai::lama::DenseMatrix<ValueType> tempDataSyn = seismogramSyntemp.getData();
        scai::lama::DenseMatrix<ValueType> tempDataObs = seismogramObstemp.getData();
        scai::lama::DenseVector<ValueType> dataTrace; 
        ValueType waterLevel = 1.0;
        Common::envelope(tempDataSyn);
        Common::envelope(tempDataObs);
        for (int i=0; i<tempDataObs.getNumRows(); i++) {
            tempDataObs.getRow(dataTrace, i);  
            if (i==0) {
                waterLevel = dataTrace.maxNorm();                
            } else {
                if (waterLevel > dataTrace.maxNorm()) waterLevel = dataTrace.maxNorm(); 
            }
        }
        waterLevel *= 1e-1;
        tempDataObs.binaryOp(tempDataSyn, scai::common::BinaryOp::SUB, tempDataObs);
        for (int i=0; i<tempDataSyn.getNumRows(); i++) {
            tempDataSyn.getRow(dataTrace, i);   
            dataTrace += waterLevel;
            tempDataSyn.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
        }
        tempDataObs.binaryOp(tempDataObs, scai::common::BinaryOp::DIVIDE, tempDataSyn);
        seismogramSyntemp.getData().binaryOp(seismogramSyntemp.getData(), scai::common::BinaryOp::MULT, tempDataObs);
        
        tempDataSyn = seismogramAdj.getData();
        Common::hilbert(tempDataSyn);
        tempDataSyn.binaryOp(tempDataSyn, scai::common::BinaryOp::MULT, tempDataObs);
        Common::hilbert(tempDataSyn);
        seismogramObstemp.getData() = tempDataSyn;  
    }
    seismogramAdj = seismogramSyntemp - seismogramObstemp;      
    
    if (seismogramSyn.getTraceType() != Acquisition::SeismogramTypeEM::HZ) {
        seismogramAdj *= -1;        
    }
}

template class KITGPI::Misfit::MisfitL2<double>;
template class KITGPI::Misfit::MisfitL2<float>;
