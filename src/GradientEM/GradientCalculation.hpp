
// #include <iosfwd>

#include <scai/common/Walltime.hpp>
#include <scai/lama.hpp>

#include <vector>

#include "GradientFactory.hpp"
#include "../ZeroLagCrossCorrelationEM/ZeroLagXcorrFactory.hpp"
#include "../Misfit/Misfit.hpp"
#include "../Preconditioning/EnergyPreconditioning.hpp"
#include "../Preconditioning/SourceReceiverTaper.hpp"
#include <AcquisitionEM/Receivers.hpp>
#include <AcquisitionEM/Sources.hpp>
#include <Configuration/Configuration.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolverEM/SourceReceiverImpl/SourceReceiverImpl.hpp>
#include <ForwardSolverEM/SourceReceiverImpl/SourceReceiverImplFactory.hpp>
#include <ForwardSolverEM/ForwardSolver.hpp>
#include <ModelparameterEM/ModelparameterFactory.hpp>
#include <WavefieldsEM/WavefieldsFactory.hpp>
#include "../WorkflowEM/Workflow.hpp"
#include "../Taper/Taper2D.hpp"
#include <Common/Hilbert.hpp>

namespace KITGPI
{
    
    /*! \brief Class to calculate the gradient for one shot
     * 
     */
    template <typename ValueType>
    class GradientCalculationEM
    {

    public:
        /* Default constructor and destructor */
        GradientCalculationEM(){};
        ~GradientCalculationEM(){};

        void allocate(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow);
        /* Calculate gradients */
        void run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolverEM<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::ReceiversEM<ValueType> &receivers, KITGPI::Acquisition::SourcesEM<ValueType> &sources, KITGPI::Acquisition::ReceiversEM<ValueType> &adjointSources, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &model, KITGPI::Gradient::GradientEM<ValueType> &gradient, std::vector<typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr> &wavefieldrecord, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, int shotNumber, KITGPI::Workflow::WorkflowEM<ValueType> const &workflow, KITGPI::Taper::Taper2D<ValueType> &wavefieldTaper2D, std::vector<typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr> &wavefieldrecordReflect);

    private:

        typedef typename KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::ZeroLagXcorrPtr ZeroLagXcorrPtr;
        ZeroLagXcorrPtr ZeroLagXcorr;

        typedef typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr wavefieldPtr;
        wavefieldPtr wavefields;
        wavefieldPtr wavefieldsTemp;

        KITGPI::Preconditioning::SourceReceiverTaper<ValueType> SourceTaper;
        KITGPI::Preconditioning::SourceReceiverTaper<ValueType> ReceiverTaper;
        KITGPI::Preconditioning::SourceReceiverTaper<ValueType> sourceReceiverTaper;
        KITGPI::Preconditioning::EnergyPreconditioning<ValueType> energyPrecond;
    };
}
