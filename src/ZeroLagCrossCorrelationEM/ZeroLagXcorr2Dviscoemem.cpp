#include "ZeroLagXcorr2Dviscoemem.hpp"

using namespace scai;

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D viscoemem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::ZeroLagXcorr2Dviscoemem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    equationType="viscoemem"; 
    numDimension=2;
    init(ctx, dist, workflow);
}

/*! \brief Initialisation which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D viscoemem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    type = equationType+std::to_string(numDimension)+"D";
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrSigmaEM, ctx, dist);
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForTauSigmaEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrEpsilonEM, ctx, dist);
    if (workflow.getInvertForEpsilonEM() || workflow.getInvertForTauEpsilonEM())
        this->initWavefield(xcorrREpsilonSigmaEM, ctx, dist);

    relaxationTime.clear();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        relaxationTime.push_back(1.0 / (2.0 * M_PI * relaxationFrequency[l]));
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::getContextPtr()
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::write(std::string filename, IndexType t, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrSigmaEM, "xcorrSigmaEM", filename, t);
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForTauSigmaEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrEpsilonEM, "xcorrEpsilonEM", filename, t);
    if (workflow.getInvertForEpsilonEM() || workflow.getInvertForTauEpsilonEM())
        this->writeWavefield(xcorrREpsilonSigmaEM, "xcorrREpsilonSigmaEM", filename, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrSigmaEM);
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForTauSigmaEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrEpsilonEM);
    if (workflow.getInvertForEpsilonEM() || workflow.getInvertForTauEpsilonEM())
        this->resetWavefield(xcorrREpsilonSigmaEM);
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
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefieldDerivative, Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //temporary wavefield allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp;    
    lama::DenseVector<ValueType> xcorrREpsilonEM;    
    lama::DenseVector<ValueType> xcorrRSigmaEM;    
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefield.getRefEX();
        xcorrSigmaEM += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefield.getRefEY();
        xcorrSigmaEM += temp;
    }
    if (workflow.getInvertForSigmaEM() || workflow.getInvertForEpsilonEM() || workflow.getInvertForTauSigmaEM() || workflow.getInvertForTauEpsilonEM() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefieldDerivative.getRefEX();
        xcorrEpsilonEM += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefieldDerivative.getRefEY();
        xcorrEpsilonEM += temp;
    }
    if (workflow.getInvertForEpsilonEM() || workflow.getInvertForTauEpsilonEM()) {
        for (int l=0; l<numRelaxationMechanisms; l++) {
            xcorrRSigmaEM = adjointWavefield.getRefRX()[l];
            xcorrRSigmaEM *= forwardWavefield.getRefRX()[l];
            temp = adjointWavefield.getRefRY()[l];
            temp *= forwardWavefield.getRefRY()[l];
            xcorrRSigmaEM += temp;
            
            xcorrREpsilonEM = adjointWavefield.getRefRX()[l];
            xcorrREpsilonEM *= forwardWavefieldDerivative.getRefRX()[l];
            temp = adjointWavefield.getRefRY()[l];
            temp *= forwardWavefieldDerivative.getRefRY()[l];
            xcorrREpsilonEM += temp; 
            
            xcorrREpsilonEM *= relaxationTime[l];
            xcorrREpsilonEM += xcorrRSigmaEM;
            xcorrREpsilonEM *= relaxationTime[l];
            xcorrREpsilonSigmaEM += xcorrREpsilonEM;
        }
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (viscoemem)
 */
template <typename ValueType>
std::string KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::getEquationType() const
{
    return (equationType);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<double>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<float>;
