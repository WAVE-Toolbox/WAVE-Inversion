#include <scai/lama.hpp>

#include "../Gradient/GradientFactory.hpp"
#include "../Workflow/Workflow.hpp"
#include "../GradientEM/GradientFactory.hpp"
#include "../WorkflowEM/Workflow.hpp"
#include "./Optimization.hpp"

namespace KITGPI
{    
    //! \brief Optimization namespace
    namespace Optimization
    {        
        /*! \brief Class to calculate the conjugate gradient direction
        * 
        */
        template <typename ValueType>
        class ConjugateGradient : public Optimization<ValueType>
        {

        public:
            /* Default constructor and destructor */
            ConjugateGradient(){};
            ~ConjugateGradient(){};

            ConjugateGradient(scai::dmemo::DistributionPtr dist);
            
            void init(scai::dmemo::DistributionPtr dist);
            void apply(KITGPI::Gradient::Gradient<ValueType> &gradient, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Configuration::Configuration config);
            void apply(KITGPI::Gradient::GradientEM<ValueType> &gradientEM, KITGPI::Workflow::WorkflowEM<ValueType> const &workflowEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelEM, KITGPI::Configuration::Configuration configEM);            

        private:
            
            void calcConjugateGradient(scai::lama::DenseVector<ValueType> &conjugateGradient, scai::lama::DenseVector<ValueType> const &gradient, scai::lama::DenseVector<ValueType> const &lastConjugateGradient, scai::lama::DenseVector<ValueType> const &lastGradient);
            
            /* Use KITGPI::Gradient::Gradient<ValueType> lastGradient and KITGPI::Gradient::Gradient<ValueType> lastConjugateGradient instead! */
            scai::lama::DenseVector<ValueType> lastGradientVp;
            scai::lama::DenseVector<ValueType> lastGradientVs;
            scai::lama::DenseVector<ValueType> lastGradientDensity;
            scai::lama::DenseVector<ValueType> lastConjugateGradientVp;
            scai::lama::DenseVector<ValueType> lastConjugateGradientVs;
            scai::lama::DenseVector<ValueType> lastConjugateGradientDensity;
            
            scai::lama::DenseVector<ValueType> lastGradientPorosity;
            scai::lama::DenseVector<ValueType> lastGradientSaturation;
            scai::lama::DenseVector<ValueType> lastGradientReflectivity;
            scai::lama::DenseVector<ValueType> lastConjugateGradientPorosity;
            scai::lama::DenseVector<ValueType> lastConjugateGradientSaturation;
            scai::lama::DenseVector<ValueType> lastConjugateGradientReflectivity;
            
            /* Use KITGPI::Gradient::GradientEM<ValueType> lastGradient and KITGPI::Gradient::GradientEM<ValueType> lastConjugateGradient instead! */
            scai::lama::DenseVector<ValueType> lastGradientSigmaEM;
            scai::lama::DenseVector<ValueType> lastGradientEpsilonEM;
            scai::lama::DenseVector<ValueType> lastGradientTauSigmaEM;
            scai::lama::DenseVector<ValueType> lastGradientTauEpsilonEM;
            scai::lama::DenseVector<ValueType> lastConjugateGradientSigmaEM;
            scai::lama::DenseVector<ValueType> lastConjugateGradientEpsilonEM;
            scai::lama::DenseVector<ValueType> lastConjugateGradientTauSigmaEM;
            scai::lama::DenseVector<ValueType> lastConjugateGradientTauEpsilonEM;
        };
    }
}
