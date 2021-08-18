#include "SteepestDescent.hpp"

/*! \brief Constructor
 * 
 * 
 \param dist
 */
template <typename ValueType>
KITGPI::Optimization::SteepestDescent<ValueType>::SteepestDescent(scai::dmemo::DistributionPtr dist)
{
    this->init(dist);   
}

/*! \brief Initialize private members with zeros
 * 
 * 
 \param dist
 */
template <typename ValueType>
void KITGPI::Optimization::SteepestDescent<ValueType>::init(scai::dmemo::DistributionPtr /*dist*/)
{

}

/*! \brief Calculate the conjugate gradient direction after Polak and Ribiere for all used parameters
 * 
 * 
 \param gradient In- and output
 \param workflow To check which parameter class is inverted for
 */
template <typename ValueType>
void KITGPI::Optimization::SteepestDescent<ValueType>::apply(KITGPI::Gradient::Gradient<ValueType> &gradient, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Configuration::Configuration config)
{    
    gradient.scale(model, workflow, config);   
}

/*! \brief Calculate the conjugate gradientEM direction after Polak and Ribiere for all used parameters
 * 
 * 
 \param gradientEM In- and output
 \param workflowEM To check which parameter class is inverted for
 */
template <typename ValueType>
void KITGPI::Optimization::SteepestDescent<ValueType>::apply(KITGPI::Gradient::GradientEM<ValueType> &gradientEM, KITGPI::Workflow::Workflow<ValueType> const &workflowEM, KITGPI::Modelparameter::Modelparameter<ValueType> const &modelEM, KITGPI::Configuration::Configuration configEM)
{    
    gradientEM.scale(modelEM, workflowEM, configEM);   
}


template class KITGPI::Optimization::SteepestDescent<double>;
template class KITGPI::Optimization::SteepestDescent<float>;
