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
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::ZeroLagXcorr2Dtmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow)
{
    equationType="tmem"; 
    numDimension=2;
    init(ctx, dist, workflow);
}


template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow)
{
    type = equationType+std::to_string(numDimension)+"D";
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrSigmaEM, ctx, dist);
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrEpsilonEM, ctx, dist);
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getContextPtr()
{
    return (xcorrEpsilonEM.getContextPtr());
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param filename file name
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        this->writeWavefield(xcorrSigmaEM, "xcorrSigmaEM", filename + type, t);
        this->writeWavefield(xcorrSigmaEMstep, "xcorrSigmaEMstep", filename + type, t);
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
    {
        this->writeWavefield(xcorrEpsilonEM, "xcorrEpsilonEM", filename + type, t);
        this->writeWavefield(xcorrEpsilonEMstep, "xcorrEpsilonEMstep", filename + type, t);
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::resetXcorr(KITGPI::Workflow::WorkflowEM<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrSigmaEM);
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrEpsilonEM);
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::update(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow, scai::IndexType gradientType, scai::IndexType decomposeType)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientType == 0 || decomposeType == 0) {   
            // Born kernel or FWI kernel
            xcorrSigmaEMstep = adjointWavefield.getRefEZ();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZ();
            xcorrSigmaEM += xcorrSigmaEMstep;           
        } else if (gradientType == 1 && decomposeType == 1) {    
            // migration kernel using up/down-going wavefields
            xcorrSigmaEMstep = adjointWavefield.getRefEZdown();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZdown();
            xcorrSigmaEM += xcorrSigmaEMstep;           
            xcorrSigmaEMstep = adjointWavefield.getRefEZup();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZup();
            xcorrSigmaEM += xcorrSigmaEMstep;      
        } else if (gradientType == 2 && decomposeType == 1) {    
            // tomographic kernel using up/down-going wavefields
            xcorrSigmaEMstep = adjointWavefield.getRefEZdown();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZup();
            xcorrSigmaEM += xcorrSigmaEMstep;           
            xcorrSigmaEMstep = adjointWavefield.getRefEZup();
            xcorrSigmaEMstep *= forwardWavefield.getRefEZdown();
            xcorrSigmaEM += xcorrSigmaEMstep;           
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        if (gradientType == 0 || decomposeType == 0) {   
            // Born kernel or FWI kernel
            xcorrEpsilonEMstep = adjointWavefield.getRefEZ();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZ();
            xcorrEpsilonEM += xcorrEpsilonEMstep;
        } else if (gradientType == 1 && decomposeType == 1) {    
            // migration kernel using up/down-going wavefields   
            xcorrEpsilonEMstep = adjointWavefield.getRefEZdown();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZdown();
            xcorrEpsilonEM += xcorrEpsilonEMstep;           
            xcorrEpsilonEMstep = adjointWavefield.getRefEZup();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZup();
            xcorrEpsilonEM += xcorrEpsilonEMstep;      
        } else if (gradientType == 2 && decomposeType == 1) {    
            // tomographic kernel using up/down-going wavefields
            xcorrEpsilonEMstep = adjointWavefield.getRefEZdown();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZup();
            xcorrEpsilonEM += xcorrEpsilonEMstep;           
            xcorrEpsilonEMstep = adjointWavefield.getRefEZup();
            xcorrEpsilonEMstep *= forwardWavefieldDerivative.getRefEZdown();
            xcorrEpsilonEM += xcorrEpsilonEMstep;      
        }
    }
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::updateHessianVectorProduct(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, Wavefields::WavefieldsEM<ValueType> &forwardWavefield2ndOrder, Wavefields::WavefieldsEM<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefield2ndOrder.getRefEZ();
        temp *= adjointWavefield.getRefEZ();
        xcorrSigmaEM += temp;
        temp = adjointWavefield2ndOrder.getRefEZ();
        temp *= forwardWavefield.getRefEZ();
        xcorrSigmaEM += temp;
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = forwardWavefield2ndOrder.getRefEZ();
        temp *= adjointWavefieldDerivative.getRefEZ();
        xcorrEpsilonEM += temp;
        temp = adjointWavefield2ndOrder.getRefEZ();
        temp *= forwardWavefieldDerivative.getRefEZ();
        xcorrEpsilonEM += temp;
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

//! \brief Getter routine for xcorrRSigmaEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getXcorrRSigmaEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrRSigmaEM in an tmem modelling")
    return (xcorrRSigmaEM);
}

//! \brief Getter routine for xcorrREpsilonEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getXcorrREpsilonEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrREpsilonEM in an tmem modelling")
    return (xcorrREpsilonEM);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<float>;
