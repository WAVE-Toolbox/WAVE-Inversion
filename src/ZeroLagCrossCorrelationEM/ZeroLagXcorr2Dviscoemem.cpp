#include "ZeroLagXcorr2Dviscoemem.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefieldsEM to zero.
 *
 * Initialisation of 2D viscoemem wavefieldsEM
 *
 \param ctx Context
 \param distEM Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::ZeroLagXcorr2Dviscoemem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    equationTypeEM="viscoemem"; 
    numDimension=2;
    init(ctx, distEM, workflowEM);
}


template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->initWavefield(xcorrSigmaEM, ctx, distEM);
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauSigmaEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->initWavefield(xcorrEpsilonEM, ctx, distEM);
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM()) {
        this->initWavefield(xcorrRSigmaEM, ctx, distEM);
        this->initWavefield(xcorrREpsilonEM, ctx, distEM);
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefieldsEM
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::getContextPtr()
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::write(std::string type, IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->writeWavefield(xcorrSigmaEM, "xcorrSigmaEM", type, t);
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauSigmaEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->writeWavefield(xcorrEpsilonEM, "xcorrEpsilonEM", type, t);
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM()) {
        this->writeWavefield(xcorrRSigmaEM, "xcorrRSigmaEM", type, t);
        this->writeWavefield(xcorrREpsilonEM, "xcorrREpsilonEM", type, t);
    }
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::writeSnapshot(IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    write(type, t, workflowEM);
}

/*! \brief Set all wavefieldsEM to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::resetXcorr(KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->resetWavefield(xcorrSigmaEM);
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauSigmaEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->resetWavefield(xcorrEpsilonEM);
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM()) {
        this->resetWavefield(xcorrRSigmaEM);
        this->resetWavefield(xcorrREpsilonEM);
    }
}

/*! \brief Function to update the result of the zero lag cross-correlation per timestep 
 * 
\begin{equation}
\label{eqn:GradientAdjoint1}
\begin{align*}
\nabla f_1(\varepsilon^\infty_e) =& \int_0^T \left( E_{1x} \pdv{E_x}{t} + E_{1y} \pdv{E_y}{t} \right) dt\\
\nabla f_1(\sigma^\infty_e) =& \int_0^T \left( E_{1x} E_x + E_{1y} E_y \right) dt\\
r_{\varepsilon e l} =& \int_0^T \left( r_{1lx} \pdv{r_{lx}}{t} + r_{1ly} \pdv{r_{ly}}{t} \right) dt\\
r_{\sigma e l} =& \int_0^T \left( r_{1lx} r_{lx} + r_{1ly} r_{ly}  \right) dt
\end{align*}
\end{equation}
 *
 *  Note that the forwardWavefieldDerivative is actually the derivative of the  forwardWavefield (see variable wavefieldrecordEM in IFOS.cpp).
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::update(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;    
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefield.getRefEX();
        xcorrSigmaEM += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefield.getRefEY();
        xcorrSigmaEM += temp;
    }
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauSigmaEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefieldDerivative.getRefEX();
        xcorrEpsilonEM += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefieldDerivative.getRefEY();
        xcorrEpsilonEM += temp;
    }
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM()) {
        temp = adjointWavefield.getRefRX();
        temp *= forwardWavefield.getRefRX();
        xcorrRSigmaEM += temp;
        temp = adjointWavefield.getRefRY();
        temp *= forwardWavefield.getRefRY();
        xcorrRSigmaEM += temp;
        
        temp = adjointWavefield.getRefRX();
        temp *= forwardWavefieldDerivative.getRefRX();
        xcorrREpsilonEM += temp;
        temp = adjointWavefield.getRefRY();
        temp *= forwardWavefieldDerivative.getRefRY();
        xcorrREpsilonEM += temp;
    }
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::updateHessianVectorProduct(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, Wavefields::WavefieldsEM<ValueType> &forwardWavefield2ndOrder, Wavefields::WavefieldsEM<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = forwardWavefield2ndOrder.getRefEX();
        temp *= adjointWavefield.getRefEX();
        xcorrSigmaEM += temp;
        temp = forwardWavefield2ndOrder.getRefEY();
        temp *= adjointWavefield.getRefEY();
        xcorrSigmaEM += temp;
        temp = adjointWavefield2ndOrder.getRefEX();
        temp *= forwardWavefield.getRefEX();
        xcorrSigmaEM += temp;
        temp = adjointWavefield2ndOrder.getRefEY();
        temp *= forwardWavefield.getRefEY();
        xcorrSigmaEM += temp;
    }
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauSigmaEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = forwardWavefield2ndOrder.getRefEX();
        temp *= adjointWavefieldDerivative.getRefEX();
        xcorrEpsilonEM += temp;
        temp = forwardWavefield2ndOrder.getRefEY();
        temp *= adjointWavefieldDerivative.getRefEY();
        xcorrEpsilonEM += temp;
        temp = adjointWavefield2ndOrder.getRefEX();
        temp *= forwardWavefieldDerivative.getRefEX();
        xcorrEpsilonEM += temp;
        temp = adjointWavefield2ndOrder.getRefEY();
        temp *= forwardWavefieldDerivative.getRefEY();
        xcorrEpsilonEM += temp;
    }
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM()) {
        temp = forwardWavefield2ndOrder.getRefRX();
        temp *= adjointWavefield.getRefRX();
        xcorrRSigmaEM += temp;
        temp = forwardWavefield2ndOrder.getRefRY();
        temp *= adjointWavefield.getRefRY();
        xcorrRSigmaEM += temp;
        temp = adjointWavefield2ndOrder.getRefRX();
        temp *= forwardWavefield.getRefRX();
        xcorrRSigmaEM += temp;
        temp = adjointWavefield2ndOrder.getRefRY();
        temp *= forwardWavefield.getRefRY();
        xcorrRSigmaEM += temp;
                
        temp = forwardWavefield2ndOrder.getRefRX();
        temp *= adjointWavefieldDerivative.getRefRX();
        xcorrREpsilonEM += temp;
        temp = forwardWavefield2ndOrder.getRefRY();
        temp *= adjointWavefieldDerivative.getRefRY();
        xcorrREpsilonEM += temp;
        temp = adjointWavefield2ndOrder.getRefRX();
        temp *= forwardWavefieldDerivative.getRefRX();
        xcorrREpsilonEM += temp;
        temp = adjointWavefield2ndOrder.getRefRY();
        temp *= forwardWavefieldDerivative.getRefRY();
        xcorrREpsilonEM += temp;
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationTypeEM (viscoemem)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::getEquationType() const
{
    return (equationTypeEM);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<float>;
