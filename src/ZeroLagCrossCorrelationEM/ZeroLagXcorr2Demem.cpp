#include "ZeroLagXcorr2Demem.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefieldsEM to zero.
 *
 * Initialisation of 2D emem wavefieldsEM
 *
 \param ctx Context
 \param distEM Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::ZeroLagXcorr2Demem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    equationTypeEM="emem"; 
    numDimension=2;
    init(ctx, distEM, workflowEM);
}


template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->initWavefield(xcorrSigmaEM, ctx, distEM);
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->initWavefield(xcorrEpsilonEM, ctx, distEM);
}

/*! \brief Returns hmemo::ContextPtr from this wavefieldsEM
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getContextPtr()
{
    return (xcorrEpsilonEM.getContextPtr());
}


/*! \brief override Methode tor write Wavefield Snapshot to file
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::write(std::string type, IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->writeWavefield(xcorrSigmaEM, "xcorrSigmaEM", type, t);
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation())
        this->writeWavefield(xcorrEpsilonEM, "xcorrEpsilonEM", type, t);
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::writeSnapshot(IndexType t, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    write(type, t, workflowEM);
}

/*! \brief Set all wavefieldsEM to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::resetXcorr(KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
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
& \nabla f_1(\varepsilon_e) = \int_0^T ( E_{1x} \pdv{E_x}{t} + E_{1y} \pdv{E_y}{t} ) dt\\
& \nabla f_1(\sigma_e) = \int_0^T ( E_{1x} E_x + E_{1y} E_y ) dt\\
\end{align*}
\end{equation}
 *
 *  Note that the forwardWavefieldDerivative is actually the derivative of the  forwardWavefield (see variable wavefieldrecordEM in IFOS.cpp).
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::update(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;    
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefield.getRefEX();
        xcorrSigmaEM += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefield.getRefEY();
        xcorrSigmaEM += temp;
    }
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefieldDerivative.getRefEX();
        xcorrEpsilonEM += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefieldDerivative.getRefEY();
        xcorrEpsilonEM += temp;
    }
}

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::updateHessianVectorProduct(Wavefields::WavefieldsEM<ValueType> &forwardWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &forwardWavefield, Wavefields::WavefieldsEM<ValueType> &adjointWavefieldDerivative, Wavefields::WavefieldsEM<ValueType> &adjointWavefield, Wavefields::WavefieldsEM<ValueType> &forwardWavefield2ndOrder, Wavefields::WavefieldsEM<ValueType> &adjointWavefield2ndOrder, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
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
    if (workflowEM.getInvertForSigmaEM() || workflowEM.getInvertForEpsilonEM() || workflowEM.getInvertForPorosity() || workflowEM.getInvertForSaturation()) {
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
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationTypeEM (emem)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getEquationType() const
{
    return (equationTypeEM);
}

//! \brief Getter routine for xcorrRSigmaEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getXcorrRSigmaEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrRSigmaEM in an emem modelling")
    return (xcorrRSigmaEM);
}

//! \brief Getter routine for xcorrREpsilonEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getXcorrREpsilonEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrREpsilonEM in an emem modelling")
    return (xcorrREpsilonEM);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<float>;
