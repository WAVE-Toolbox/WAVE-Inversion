#include <scai/lama.hpp>

#include "../Gradient/GradientFactory.hpp"
#include "../Workflow/Workflow.hpp"
#include "./Optimization.hpp"

namespace KITGPI
{    
    //! \brief Optimization namespace
    namespace Optimization
    {        
        /*! \brief Class to calculate the conjugate gradientEM direction
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
            scai::lama::DenseVector<ValueType> lastConjugateGradientPorosity;
            scai::lama::DenseVector<ValueType> lastConjugateGradientSaturation;
        };
    }
}
