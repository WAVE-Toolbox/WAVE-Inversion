#include "ZeroLagXcorr2Dviscotmem.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefieldsEM to zero.
 *
 * Initialisation of 2D viscotmem wavefieldsEM
 *
 \param ctx Context
 \param distEM Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<ValueType>::ZeroLagXcorr2Dviscotmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    equationTypeEM="viscotmem"; 
    numDimension=2;
    init(ctx, distEM, workflowEM);
}


template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
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
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<ValueType>::getContextPtr()
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<ValueType>::write(std::string type, IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<ValueType>::writeSnapshot(IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    write(type, t, workflowEM);
}

/*! \brief Set all wavefieldsEM to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<ValueType>::resetXcorr(KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
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
\nabla f_1(\varepsilon^\infty_e) =& \int_0^T \left( E_{1z} \pdv{E_z}{t} \right) dt\\
\nabla f_1(\sigma^\infty_e) =& \int_0^T \left( E_{1z} E_z  \right) dt\\
r_{\varepsilon e l} =& \int_0^T \left( r_{1lz} \pdv{r_{lz}}{t} \right) dt\\
r_{\sigma e l} =& \int_0^T \left( r_{1lz} r_{lz}  \right) dt
\end{align*}
\end{equation}
 *
 *  Note that the forwardWavefieldDerivative is actually the derivative of the  forwardWavefield (see variable wavefieldrecordEM in IFOS.cpp).
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<ValueType>::update(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;    
//     std::cout << "update adjointWavefield.getRefEZ() = " << adjointWavefield.getRefEZ() << "\n" << std::endl;
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEZ();
        temp *= forwardWavefield.getRefEZ();
        xcorrSigmaEM += temp;
    }
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauSigmaEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEZ();
        temp *= forwardWavefieldDerivative.getRefEZ();
        xcorrEpsilonEM += temp;
    }
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM()) {
        temp = adjointWavefield.getRefRZ();
        temp *= forwardWavefield.getRefRZ();
        xcorrRSigmaEM += temp;
        
        temp = adjointWavefield.getRefRZ();
        temp *= forwardWavefieldDerivative.getRefRZ();
        xcorrREpsilonEM += temp;
    }
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<ValueType>::updateHessianVectorProduct(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, Wavefields::WavefieldsEM<ValueType> &forwardWavefield2ndOrder, Wavefields::WavefieldsEM<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = forwardWavefield2ndOrder.getRefEZ();
        temp *= adjointWavefield.getRefEZ();
        xcorrSigmaEM += temp;
        temp = adjointWavefield2ndOrder.getRefEZ();
        temp *= forwardWavefield.getRefEZ();
        xcorrSigmaEM += temp;
    }
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauSigmaEM() || workflowEM.getInvertForTauEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = forwardWavefield2ndOrder.getRefEZ();
        temp *= adjointWavefieldDerivative.getRefEZ();
        xcorrEpsilonEM += temp;
        temp = adjointWavefield2ndOrder.getRefEZ();
        temp *= forwardWavefieldDerivative.getRefEZ();
        xcorrEpsilonEM += temp;
    }
    if (workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForTauEpsilonEM()) {
        temp = forwardWavefield2ndOrder.getRefRZ();
        temp *= adjointWavefield.getRefRZ();
        xcorrRSigmaEM += temp;
        temp = adjointWavefield2ndOrder.getRefRZ();
        temp *= forwardWavefield.getRefRZ();
        xcorrRSigmaEM += temp;
                
        temp = forwardWavefield2ndOrder.getRefRZ();
        temp *= adjointWavefieldDerivative.getRefRZ();
        xcorrREpsilonEM += temp;
        temp = adjointWavefield2ndOrder.getRefRZ();
        temp *= forwardWavefieldDerivative.getRefRZ();
        xcorrREpsilonEM += temp;
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationTypeEM (viscotmem)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<ValueType>::getEquationType() const
{
    return (equationTypeEM);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscotmem<float>;
