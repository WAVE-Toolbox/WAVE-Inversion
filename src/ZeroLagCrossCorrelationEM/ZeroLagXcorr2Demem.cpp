#include "ZeroLagXcorr2Demem.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D emem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::ZeroLagXcorr2Demem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
{
    equationType="emem"; 
    numDimension=2;
    init(ctx, dist, workflow, config);
}


template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config)
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
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getContextPtr()
{
    return (xcorrEpsilonEM.getContextPtr());
}


/*! \brief override Method to write Wavefield Snapshot to file
 \param filename file name
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrSigmaEM, "xcorrSigmaEM", filename, t);
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrEpsilonEM, "xcorrEpsilonEM", filename, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
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
& \nabla f_1(\varepsilon_e) = \int_0^T ( E_{1x} \pdv{E_x}{t} + E_{1y} \pdv{E_y}{t} ) dt\\
& \nabla f_1(\sigma_e) = \int_0^T ( E_{1x} E_x + E_{1y} E_y ) dt\\
\end{align*}
\end{equation}
 *
 *  Note that the forwardWavefieldDerivative is actually the derivative of the  forwardWavefield (see variable wavefieldrecordEM in IFOS.cpp).
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;    
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefield.getRefEX();
        xcorrSigmaEM += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefield.getRefEY();
        xcorrSigmaEM += temp;
        if (gradientKernel != 0) {           
        }
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefieldDerivative.getRefEX();
        xcorrEpsilonEM += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefieldDerivative.getRefEY();
        xcorrEpsilonEM += temp;
        if (gradientKernel != 0) {     
        }
    }
}

/*! \brief Gather wavefields in the time domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::gatherWavefields(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep)
{
}

/*! \brief Sum wavefields in the frequency domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::sumWavefields(KITGPI::Workflow::Workflow<ValueType> const &workflow, ValueType DT)
{
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (emem)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Getter routine for xcorrREpsilonSigmaEM
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<ValueType>::getXcorrREpsilonSigmaEM() const
{
    COMMON_THROWEXCEPTION("There is no xcorrREpsilonSigmaEM in an emem modelling")
    return (xcorrREpsilonSigmaEM);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Demem<float>;
