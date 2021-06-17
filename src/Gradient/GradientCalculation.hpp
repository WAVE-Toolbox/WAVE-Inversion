
// #include <iosfwd>

#include <scai/common/Walltime.hpp>
#include <scai/lama.hpp>

#include <vector>

#include "GradientFactory.hpp"
#include "../ZeroLagCrossCorrelation/ZeroLagXcorrFactory.hpp"
#include "../Misfit/Misfit.hpp"
#include "../Preconditioning/EnergyPreconditioning.hpp"
#include "../Preconditioning/SourceReceiverTaper.hpp"
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <Configuration/Configuration.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>
#include "../Workflow/Workflow.hpp"
#include "../Taper/Taper2D.hpp"

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

        void allocate(KITGPI::Configuration::Configuration config, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::Workflow::Workflow<ValueType> const &workflow);
        /* Calculate gradients */
        void run(scai::dmemo::CommunicatorPtr commAll, KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> const &adjointSources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Gradient::Gradient<ValueType> &gradient, std::vector<typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr> &wavefieldrecord, KITGPI::Configuration::Configuration config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, int shotNumber, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Taper::Taper2D<ValueType> wavefieldTaper2D);

    private:

        typedef typename KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType>::ZeroLagXcorrPtr ZeroLagXcorrPtr;
        ZeroLagXcorrPtr ZeroLagXcorr;

        typedef typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldPtr;
        wavefieldPtr wavefields;
        wavefieldPtr wavefieldsTemp;

        KITGPI::Preconditioning::SourceReceiverTaper<ValueType> SourceTaper;
        KITGPI::Preconditioning::SourceReceiverTaper<ValueType> ReceiverTaper;
        KITGPI::Preconditioning::EnergyPreconditioning<ValueType> energyPrecond;
    };
}
