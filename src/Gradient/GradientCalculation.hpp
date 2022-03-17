
// #include <iosfwd>

#include <scai/common/Walltime.hpp>
#include <scai/lama.hpp>

#include <vector>

#include "GradientFactory.hpp"
#include "../ZeroLagCrossCorrelation/ZeroLagXcorrFactory.hpp"
#include "../Misfit/Misfit.hpp"
#include "../Misfit/MisfitFactory.hpp"
#include "../Preconditioning/EnergyPreconditioning.hpp"
#include "../Preconditioning/SourceReceiverTaper.hpp"
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <Configuration/Configuration.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolver/SourceReceiverImpl/SourceReceiverImplFactory.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>
#include "../Workflow/Workflow.hpp"
#include "../Taper/Taper2D.hpp"

using namespace scai;

namespace KITGPI
{    
    /*! \brief Class to calculate the gradient for one shot
     * 
     */
    template <typename ValueType>
    class GradientCalculation
    {

    public:
        /* Default constructor and destructor */
        GradientCalculation(){};
        ~GradientCalculation(){};

        void allocate(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distInversion, scai::hmemo::ContextPtr ctx, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType numShotPerSuperShot);
        void gatherWavefields(KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsInversion, scai::lama::DenseVector<ValueType> sourceFC, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::IndexType tStep, ValueType DT, bool isAdjoint = false, bool isReflect = false);
        
        /* Calculate gradients */
        void run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> sources, KITGPI::Acquisition::Receivers<ValueType> const &adjointSources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, std::vector<typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr> &wavefieldrecord, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, int shotNumber, int shotIndTrue, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Taper::Taper2D<ValueType> wavefieldTaper2D, std::vector<typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr> &wavefieldrecordReflect, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, KITGPI::Preconditioning::EnergyPreconditioning<ValueType> energyPrecond, KITGPI::Preconditioning::EnergyPreconditioning<ValueType> energyPrecondReflect, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode);

    private:

        typedef typename KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::ZeroLagXcorrPtr ZeroLagXcorrPtr;
        ZeroLagXcorrPtr ZeroLagXcorr;
        ZeroLagXcorrPtr ZeroLagXcorrReflect;

        typedef typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
        wavefieldPtr wavefields;
        wavefieldPtr wavefieldsReflect;
        wavefieldPtr wavefieldsTemp;
        wavefieldPtr wavefieldsAdjointTemp;

        KITGPI::Preconditioning::SourceReceiverTaper<ValueType> SourceTaper;
        KITGPI::Preconditioning::SourceReceiverTaper<ValueType> ReceiverTaper;
        KITGPI::Preconditioning::SourceReceiverTaper<ValueType> sourceReceiverTaper;
    };
}
