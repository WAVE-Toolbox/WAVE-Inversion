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
        class SteepestDescent : public Optimization<ValueType>
        {

        public:
            /* Default constructor and destructor */
            SteepestDescent(){};
            ~SteepestDescent(){};

            SteepestDescent(scai::dmemo::DistributionPtr distEM);
            
            void init(scai::dmemo::DistributionPtr dist);
            void apply(KITGPI::Gradient::Gradient<ValueType> &gradient, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Configuration::Configuration config);

        };
    }
}
