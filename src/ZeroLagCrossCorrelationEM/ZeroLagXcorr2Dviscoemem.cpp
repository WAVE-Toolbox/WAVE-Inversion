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
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::ZeroLagXcorr2Dviscoemem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot)
{
    equationType="viscoemem"; 
    numDimension=2;
    init(ctx, dist, workflow, config, numShotPerSuperShot);
}

/*! \brief Initialisation which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D viscoemem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Configuration::Configuration config, scai::IndexType numShotPerSuperShot)
{
    type = equationType+std::to_string(numDimension)+"D";
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrSigma, ctx, dist);
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrEpsilon, ctx, dist);
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->initWavefield(xcorrREpsilonSigma, ctx, dist);

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
    return (xcorrEpsilon.getContextPtr());
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
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrSigma, "xcorrSigma", filename, t);
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrEpsilon, "xcorrEpsilon", filename, t);
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->writeWavefield(xcorrREpsilonSigma, "xcorrREpsilonSigma", filename, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::resetXcorr(KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrSigma);
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrEpsilon);
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation())
        this->resetWavefield(xcorrREpsilonSigma);
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
    lama::DenseVector<ValueType> xcorrRSigma;    
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefield.getRefEX();
        xcorrSigma += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefield.getRefEY();
        xcorrSigma += temp;
    }
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        temp = adjointWavefield.getRefEX();
        temp *= forwardWavefieldDerivative.getRefEX();
        xcorrEpsilon += temp;
        temp = adjointWavefield.getRefEY();
        temp *= forwardWavefieldDerivative.getRefEY();
        xcorrEpsilon += temp;
    }
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        for (int l=0; l<numRelaxationMechanisms; l++) {
            xcorrRSigma = adjointWavefield.getRefRX()[l];
            xcorrRSigma *= forwardWavefield.getRefEX();
            temp = adjointWavefield.getRefRY()[l];
            temp *= forwardWavefield.getRefEY();
            xcorrRSigma += temp; 
            
            xcorrREpsilonSigma += xcorrRSigma;
        }
    }
}

/*! \brief Gather wavefields in the time domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::gatherWavefields(Wavefields::Wavefields<ValueType> &wavefields, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint)
{
}

/*! \brief Sum wavefields in the frequency domain
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::sumWavefields(scai::dmemo::CommunicatorPtr commShot, std::string filename, IndexType snapType, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::lama::DenseVector<ValueType> sourceFC, ValueType DT, scai::IndexType shotNumber, std::vector<scai::lama::SparseVector<ValueType>> taperEncode)
{
}

/*! \brief Recover ZeroLagXcorr to original size
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dviscoemem<ValueType>::applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        xcorrSigma = lhs * xcorrSigma;
    }
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        xcorrEpsilon = lhs * xcorrEpsilon;
    }
    if (workflow.getInvertForSigma() || workflow.getInvertForEpsilon() || workflow.getInvertForTauSigma() || workflow.getInvertForTauEpsilon() || workflow.getInvertForPorosity() || workflow.getInvertForSaturation()) {
        xcorrREpsilonSigma = lhs * xcorrREpsilonSigma;
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
