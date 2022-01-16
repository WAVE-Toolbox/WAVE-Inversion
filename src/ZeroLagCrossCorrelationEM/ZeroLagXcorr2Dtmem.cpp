#include "ZeroLagXcorr2Dtmem.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D tmem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::ZeroLagXcorr2Dtmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    equationType="tmem"; 
    numDimension=2;
    init(ctx, dist, workflow);
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    type = equationType+std::to_string(numDimension)+"D";
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrSigmaEM, ctx, dist);
        if (gradientType != 0 && decomposition != 0) {
            this->initWavefield(xcorrSigmaEMSuRu, ctx, dist);
            this->initWavefield(xcorrSigmaEMSdRd, ctx, dist);
            this->initWavefield(xcorrSigmaEMSuRd, ctx, dist);
            this->initWavefield(xcorrSigmaEMSdRu, ctx, dist);
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->initWavefield(xcorrEpsilonEM, ctx, dist);
        if (gradientType != 0 && decomposition != 0) {
            this->initWavefield(xcorrEpsilonEMSuRu, ctx, dist);
            this->initWavefield(xcorrEpsilonEMSdRd, ctx, dist);
            this->initWavefield(xcorrEpsilonEMSuRd, ctx, dist);
            this->initWavefield(xcorrEpsilonEMSdRu, ctx, dist);
        }
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getContextPtr()
{
    return (xcorrEpsilonEM.getContextPtr());
}

/*! \brief override Method to write Wavefield Snapshot to file
 *
 *
 \param filename file name
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrSigmaEM, "xcorrSigmaEM", filename, t);
        if (decomposition == 0) {
            this->writeWavefield(xcorrSigmaEMstep, "xcorrSigmaEM.step", filename, t);
        } else if (gradientType == 1 && decomposition == 1) {
            this->writeWavefield(xcorrSigmaEMSdRu, "xcorrSigmaEM.SdRu", filename, t);
            this->writeWavefield(xcorrSigmaEMSuRd, "xcorrSigmaEM.SuRd", filename, t);
        } else if (gradientType == 2 && decomposition == 1) {
            this->writeWavefield(xcorrSigmaEMSuRu, "xcorrSigmaEM.SuRu", filename, t);
            this->writeWavefield(xcorrSigmaEMSdRd, "xcorrSigmaEM.SdRd", filename, t);
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
    {
        this->writeWavefield(xcorrEpsilonEM, "xcorrEpsilonEM", filename, t);
        if (decomposition == 0) {
            this->writeWavefield(xcorrEpsilonEMstep, "xcorrEpsilonEM.step", filename, t);
        } else if (gradientType == 1 && decomposition == 1) {
            this->writeWavefield(xcorrEpsilonEMSdRu, "xcorrEpsilonEM.SdRu", filename, t);
            this->writeWavefield(xcorrEpsilonEMSuRd, "xcorrEpsilonEM.SuRd", filename, t);
        } else if (gradientType == 2 && decomposition == 1) {
            this->writeWavefield(xcorrEpsilonEMSuRu, "xcorrEpsilonEM.SuRu", filename, t);
            this->writeWavefield(xcorrEpsilonEMSdRd, "xcorrEpsilonEM.SdRd", filename, t);
        }
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrSigmaEM);
        if (gradientType != 0 && decomposition != 0) {
            this->resetWavefield(xcorrSigmaEMSuRu);
            this->resetWavefield(xcorrSigmaEMSdRd);
            this->resetWavefield(xcorrSigmaEMSuRd);
            this->resetWavefield(xcorrSigmaEMSdRu);
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->resetWavefield(xcorrEpsilonEM);
        if (gradientType != 0 && decomposition != 0) {
            this->resetWavefield(xcorrEpsilonEMSuRu);
            this->resetWavefield(xcorrEpsilonEMSdRd);
            this->resetWavefield(xcorrEpsilonEMSuRd);
            this->resetWavefield(xcorrEpsilonEMSdRu);
        }
    }
}

/*! \brief Function to update the result of the zero lag cross-correlation per timestep 
 * 
\begin{equation}
\label{eqn:GradientAdjoint1}
\begin{align*}
& \nabla f_1(\varepsilon_e) = \int_0^T ( E_{1z} \pdv{E_z}{t} ) dt\\
& \nabla f_1(\sigma_e) = \int_0^T ( E_{1z} E_z  ) dt\\
\end{align*}
\end{equation}
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientType == 0 || decomposition == 0) {   
            // Born kernel or FWI kernel
            xcorrSigmaEMstep = adjointWavefield.getRefEZ();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZ();
            xcorrSigmaEM += xcorrSigmaEMstep;           
        } else if (gradientType == 1 && decomposition == 1) {    
            // migration kernel using up/down-going wavefields
            xcorrSigmaEMstep = adjointWavefield.getRefEZdown();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZdown();
            xcorrSigmaEMSdRu += xcorrSigmaEMstep;           
            xcorrSigmaEMstep = adjointWavefield.getRefEZup();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZup();
            xcorrSigmaEMSuRd += xcorrSigmaEMstep;      
            xcorrSigmaEM = xcorrSigmaEMSdRu + xcorrSigmaEMSuRd;    
        } else if (gradientType == 2 && decomposition == 1) {    
            // tomographic kernel using up/down-going wavefields
            xcorrSigmaEMstep = adjointWavefield.getRefEZdown();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZup();
            xcorrSigmaEMSuRu += xcorrSigmaEMstep;           
            xcorrSigmaEMstep = adjointWavefield.getRefEZup();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZdown();
            xcorrSigmaEMSdRd += xcorrSigmaEMstep;  
            xcorrSigmaEM = xcorrSigmaEMSuRu + xcorrSigmaEMSdRd;           
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientType == 0 || decomposition == 0) {   
            // Born kernel or FWI kernel
            xcorrEpsilonEMstep = adjointWavefield.getRefEZ();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZ();
            xcorrEpsilonEM += xcorrEpsilonEMstep;
        } else if (gradientType == 1 && decomposition == 1) {    
            // migration kernel using up/down-going wavefields   
            xcorrEpsilonEMstep = adjointWavefield.getRefEZdown();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZdown();
            xcorrEpsilonEMSdRu += xcorrEpsilonEMstep;           
            xcorrEpsilonEMstep = adjointWavefield.getRefEZup();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZup();
            xcorrEpsilonEMSuRd += xcorrEpsilonEMstep;      
            xcorrEpsilonEM = xcorrEpsilonEMSdRu + xcorrEpsilonEMSuRd;   
        } else if (gradientType == 2 && decomposition == 1) {    
            // tomographic kernel using up/down-going wavefields
            xcorrEpsilonEMstep = adjointWavefield.getRefEZdown();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZup();
            xcorrEpsilonEMSuRu += xcorrEpsilonEMstep;           
            xcorrEpsilonEMstep = adjointWavefield.getRefEZup();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZdown();
            xcorrEpsilonEMSdRd += xcorrEpsilonEMstep;  
            xcorrEpsilonEM = xcorrEpsilonEMSuRu + xcorrEpsilonEMSdRd;      
        }
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (tmem)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Getter routine for xcorrREpsilonSigmaEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getXcorrREpsilonSigmaEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrREpsilonSigmaEM in an tmem modelling")
    return (xcorrREpsilonSigmaEM);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<float>;
