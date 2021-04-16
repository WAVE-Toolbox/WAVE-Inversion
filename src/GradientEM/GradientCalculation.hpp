
// #include <iosfwd>

#include <scai/common/Walltime.hpp>
#include <scai/lama.hpp>

#include <vector>

#include "GradientFactory.hpp"
#include "../ZeroLagCrossCorrelationEM/ZeroLagXcorrFactory.hpp"
#include "../Misfit/Misfit.hpp"
#include "../Preconditioning/SourceReceiverTaper.hpp"
#include <AcquisitionEM/Receivers.hpp>
#include <AcquisitionEM/Sources.hpp>
#include <Configuration/Configuration.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolverEM/SourceReceiverImpl/SourceReceiverImpl.hpp>
#include <ForwardSolverEM/SourceReceiverImpl/FDTD2Demem.hpp>
#include <ForwardSolverEM/ForwardSolver.hpp>
#include <ModelparameterEM/ModelparameterFactory.hpp>
#include <WavefieldsEM/WavefieldsFactory.hpp>
#include "../WorkflowEM/Workflow.hpp"
#include "../Taper/Taper2D.hpp"

namespace KITGPI
{
    
    /*! \brief Class to calculate the gradientEM for one shot
     * 
     */
    template <typename ValueType>
    class GradientCalculationEM
    {

    public:
        /* Default constructor and destructor */
        GradientCalculationEM(){};
        ~GradientCalculationEM(){};

        void allocate(KITGPI::Configuration::Configuration configEM, scai::dmemo::DistributionPtr distEM, scai::hmemo::ContextPtr ctx, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM);
        /* Calculate gradients */
        void run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolverEM<ValueType> &solverEM, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivativesEM, KITGPI::Acquisition::ReceiversEM<ValueType> &ReceiversEM, KITGPI::Acquisition::SourcesEM<ValueType> &sourcesEM, KITGPI::Acquisition::ReceiversEM<ValueType> &adjointSourcesEM, KITGPI::Modelparameter::ModelparameterEM<ValueType> const &modelEM, KITGPI::Gradient::GradientEM<ValueType> &gradientEM, std::vector<typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr> &wavefieldrecordEM, KITGPI::Configuration::Configuration configEM, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinatesEM, int shotNumber, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM, KITGPI::Taper::Taper2D<ValueType> &wavefieldTaper2DEM);

    private:

        typedef typename KITGPI::ZeroLagXcorr::ZeroLagXcorrEM<ValueType>::ZeroLagXcorrPtr ZeroLagXcorrPtr;
        ZeroLagXcorrPtr ZeroLagXcorr;

        typedef typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr wavefieldPtrEM;
        wavefieldPtrEM wavefieldsEM;
        wavefieldPtrEM wavefieldsTempEM;

        KITGPI::Preconditioning::SourceReceiverTaper <ValueType> SourceTaper;

    };
}
