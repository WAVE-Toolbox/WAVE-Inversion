#include <scai/lama.hpp>

#include "../Gradient/GradientFactory.hpp"
#include "../Workflow/Workflow.hpp"

namespace KITGPI
{
    
    /*! \brief Class to calculate the conjugate gradient direction
     * 
     */
    template <typename ValueType>
    class ConjugateGradient
    {

    public:
        /* Default constructor and destructor */
        ConjugateGradient(){};
        ~ConjugateGradient(){};

        ConjugateGradient(scai::dmemo::DistributionPtr dist);
        
        void init(scai::dmemo::DistributionPtr dist);
        void calc(KITGPI::Gradient::Gradient<ValueType> &gradient, KITGPI::Workflow::Workflow<ValueType> const &workflow, int iteration);
        

    private:
        
        void calcConjugateGradient(scai::lama::DenseVector<ValueType> &conjugateGradient, scai::lama::DenseVector<ValueType> const &gradient, scai::lama::DenseVector<ValueType> const &lastConjugateGradient, scai::lama::DenseVector<ValueType> const &lastGradient);
        
        /* Use KITGPI::Gradient::Gradient<ValueType> lastGradient and KITGPI::Gradient::Gradient<ValueType> lastConjugateGradient instead! */
        scai::lama::DenseVector<ValueType> lastGradientVp;
        scai::lama::DenseVector<ValueType> lastGradientVs;
        scai::lama::DenseVector<ValueType> lastGradientDensity;
        scai::lama::DenseVector<ValueType> lastConjugateGradientVp;
        scai::lama::DenseVector<ValueType> lastConjugateGradientVs;
        scai::lama::DenseVector<ValueType> lastConjugateGradientDensity;

    };
}
