#pragma once

#include <scai/lama.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>

#include "../Gradient/GradientFactory.hpp"
#include "../Workflow/Workflow.hpp"


namespace KITGPI
{

    //! \brief Optimization namespace
    namespace Optimization
    {

        /*! \brief Abstract class to optimize the gradient
         *
	     * This class implements some methods which are required by all derived classes.
         * As this class is an abstract class, all constructors are protected.
         */
        template <typename ValueType>
        class Optimization
        {
            
          public:
              
              //! \brief Optimization pointer
              typedef std::shared_ptr<Optimization<ValueType>> OptimizationPtr;
            
              virtual void init(scai::dmemo::DistributionPtr dist) = 0;
              virtual void apply(KITGPI::Gradient::Gradient<ValueType> &gradient, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Modelparameter::Modelparameter<ValueType> const &model) = 0;
            
	    
          protected:
              
              //! Default constructor and destructor.
              Optimization(){};
              ~Optimization(){};
                            
        };
    }
}
