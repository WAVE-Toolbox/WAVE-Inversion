#include "ZeroLagXcorr2Delastic.hpp"

using namespace scai;

template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    invertForVp = config.get<bool>("invertForVp");
    invertForVs = config.get<bool>("invertForVs");
    invertForDensity = config.get<bool>("invertForDensity");
    
    if (invertForDensity)
        this->initWavefield(VSum, ctx, dist);

    if ((invertForVp) || (invertForVs) || (invertForDensity)) 
        this->initWavefield(NormalStressSum, ctx, dist);

    if ((invertForVs) || (invertForDensity)){
        this->initWavefield(ShearStress, ctx, dist);
        this->initWavefield(NormalStressDiff, ctx, dist);
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::getContextPtr()
{
    return (VSum.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D elastic wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::ZeroLagXcorr2Delastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(config, ctx, dist);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::write(std::string type, IndexType t)
{
    this->writeWavefield(VSum, "VSum", type, t);
    COMMON_THROWEXCEPTION("write correlated elastic wavefields is not implemented yet.")
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::writeSnapshot(IndexType t)
{
    write(type, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::reset()
{
    if (invertForDensity)
    this->resetWavefield(VSum);

    if ((invertForVs)|| (invertForDensity)) {
    this->resetWavefield(ShearStress);
    this->resetWavefield(NormalStressDiff);
    }
    
    if ((invertForVp) || (invertForVs)|| (invertForDensity)){
    this->resetWavefield(NormalStressSum);
    }
}

/*! \brief function to update the result of the zero lag cross-correlation for per timestep 
 */
template <typename ValueType>
void KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::update(Wavefields::Wavefields<ValueType> &forwardWavefield, Wavefields::Wavefields<ValueType> &adjointWavefield)
{
    //temporary wavefields allocated for every timestep (might be inefficient)
    lama::DenseVector<ValueType> temp1;
    lama::DenseVector<ValueType> temp2;
    if ((invertForVp) || (invertForVs) || (invertForDensity)) {
        temp1 = forwardWavefield.getRefSxx() + forwardWavefield.getRefSyy();
        temp2 = adjointWavefield.getRefSxx() + adjointWavefield.getRefSyy();
        temp1 *= temp2;
        NormalStressSum += temp1;
    }

    if ((invertForVs) || (invertForDensity)){
        temp1 = forwardWavefield.getRefSxx() - forwardWavefield.getRefSyy();
        temp2 = adjointWavefield.getRefSxx() - adjointWavefield.getRefSyy();
        temp1 *= temp2;
        NormalStressDiff += temp1;
        temp1 = forwardWavefield.getRefSxy();
        temp1 *= adjointWavefield.getRefSxy();
        ShearStress += temp1;
    }

    if (invertForDensity) {
        temp1 = forwardWavefield.getRefVX();
        temp1 *= adjointWavefield.getRefVX();
        VSum += temp1;
        temp1 = forwardWavefield.getRefVY();
        temp1 *= adjointWavefield.getRefVY();
        VSum += temp1;
    }
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::getP() const
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D elastic case.");
    return (P);
}

template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<float>;
template class KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<double>;
