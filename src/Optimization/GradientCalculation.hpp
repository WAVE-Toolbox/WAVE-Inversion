
// #include <iosfwd>

#include <scai/common/Walltime.hpp>
#include <scai/lama.hpp>

#include <vector>

#include "../Gradient/GradientFactory.hpp"
#include "./ZeroLagCrossCorrelation/ZeroLagXcorrFactory.hpp"
#include "Misfit/Misfit.hpp"
#include "SourceReceiverTaper.hpp"
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <Configuration/Configuration.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>


/*! \brief The class GradientCalculation holds the gradients for inversion
 *
 */
template <typename ValueType>
class GradientCalculation
{

  public:
    /* Default constructor and destructor */
    GradientCalculation(){};
    ~GradientCalculation(){};

    void allocate(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx);
    /* Calculate gradients */
    void run(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> const &adjointSources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Gradient::Gradient<ValueType> &gradient, std::vector<typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr> &wavefieldrecord, KITGPI::Configuration::Configuration config, IndexType iteration, int shotNumber);

  private:

    typedef typename KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::ZeroLagXcorrPtr ZeroLagXcorrPtr;
    ZeroLagXcorrPtr ZeroLagXcorr;

    typedef typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
    wavefieldPtr wavefields;

    SourceReceiverTaper <ValueType> SourceTaper;

};
