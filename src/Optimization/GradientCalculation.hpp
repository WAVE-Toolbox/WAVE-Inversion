
// #include <iosfwd>

#include <scai/common/Walltime.hpp>
#include <scai/lama.hpp>

#include "../Gradient/GradientFactory.hpp"
#include "./ZeroLagCrossCorrelation/ZeroLagXcorrFactory.hpp"
#include "Misfit.hpp"
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <Configuration/Configuration.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>

template <typename ValueType>
class GradientCalculation
{

  public:
    /* Default constructor and destructor */
    GradientCalculation(){};
    ~GradientCalculation(){};

    void allocate(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr no_dist_NT, scai::dmemo::CommunicatorPtr comm, scai::hmemo::ContextPtr ctx);
    /* Calculate gradients */
    void calc(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Gradient::Gradient<ValueType> &gradient, KITGPI::Configuration::Configuration config, IndexType iteration, Misfit<ValueType> &dataMisfit);

  private:
    typedef typename KITGPI::Gradient::Gradient<ValueType>::GradientPtr GradientPtr;
    GradientPtr GradientPerShot;

    typedef typename KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::ZeroLagXcorrPtr ZeroLagXcorrPtr;
    ZeroLagXcorrPtr ZeroLagXcorr;

    typedef typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
    wavefieldPtr wavefields;

    std::vector<wavefieldPtr> wavefieldrecord;
};
