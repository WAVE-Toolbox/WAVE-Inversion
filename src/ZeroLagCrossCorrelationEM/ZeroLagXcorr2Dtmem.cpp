#include "ZeroLagXcorr2Dtmem.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefieldsEM to zero.
 *
 * Initialisation of 2D tmem wavefieldsEM
 *
 \param ctx Context
 \param distEM Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::ZeroLagXcorr2Dtmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    equationTypeEM="tmem"; 
    numDimension=2;
    init(ctx, distEM, workflowEM);
}


template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->initWavefield(xcorrSigmaEM, ctx, distEM);
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->initWavefield(xcorrEpsilonEM, ctx, distEM);
}

/*! \brief Returns hmemo::ContextPtr from this wavefieldsEM
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getContextPtr()
{
    return (xcorrEpsilonEM.getContextPtr());
}


/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::write(std::string type, IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->writeWavefield(xcorrSigmaEM, "xcorrSigmaEM", type, t);
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->writeWavefield(xcorrEpsilonEM, "xcorrEpsilonEM", type, t);
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 * \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::writeSnapshot(IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    write(type, t, workflowEM);
}

/*! \brief Set all wavefieldsEM to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::resetXcorr(KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->resetWavefield(xcorrSigmaEM);
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
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
 *
 *  Note that the forwardWavefieldDerivative is actually the derivative of the  forwardWavefield (see variable wavefieldrecordEM in IFOS.cpp).
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::update(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;    
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEZ();
        temp *= forwardWavefield.getRefEZ();
        xcorrSigmaEM += temp;
    }
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEZ();
        temp *= forwardWavefieldDerivative.getRefEZ();
        xcorrEpsilonEM += temp;
    }
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::updateHessianVectorProduct(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, Wavefields::WavefieldsEM<ValueType> &forwardWavefield2ndOrder, Wavefields::WavefieldsEM<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = forwardWavefield2ndOrder.getRefEZ();
        temp *= adjointWavefield.getRefEZ();
        xcorrSigmaEM += temp;
        temp = adjointWavefield2ndOrder.getRefEZ();
        temp *= forwardWavefield.getRefEZ();
        xcorrSigmaEM += temp;
    }
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
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

/*! \brief Get equationTypeEM (tmem)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dtmem<ValueType>::getEquationType() const
{
    return (equationTypeEM);
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
