#include "ConjugateGradient.hpp"

/*! \brief Constructor
 * 
 * 
 \param dist
 */
template <typename ValueType>
KITGPI::Optimization::ConjugateGradient<ValueType>::ConjugateGradient(scai::dmemo::DistributionPtr dist)
{
    this->init(dist);
   
}

/*! \brief Initialize private members with zeros
 * 
 * 
 \param dist
 */
template <typename ValueType>
void KITGPI::Optimization::ConjugateGradient<ValueType>::init(scai::dmemo::DistributionPtr dist)
{
    lastGradientVp.setSameValue(dist, 0 );
    lastGradientVs.setSameValue(dist, 0 );
    lastGradientDensity.setSameValue(dist, 0 );
    lastConjugateGradientVp.setSameValue(dist, 0 );
    lastConjugateGradientVs.setSameValue(dist, 0 );
    lastConjugateGradientDensity.setSameValue(dist, 0 );
   
}


/*! \brief Calculate the conjugate gradient direction after Polak and Ribiere for all used parameters
 * 
 * 
 \param gradient In- and output
 \param workflow To check which parameter class is inverted for
 */
template <typename ValueType>
void KITGPI::Optimization::ConjugateGradient<ValueType>::apply(KITGPI::Gradient::Gradient<ValueType> &gradient, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Modelparameter::Modelparameter<ValueType> const &model)
{
    
    /* Should an automatic direction reset be implemented? -> beta = max{beta^PR, 0} */
    
    if(workflow.iteration==0){
        
        if(workflow.getInvertForVp() == 1){
            lastGradientVp = gradient.getVelocityP();
            lastConjugateGradientVp = gradient.getVelocityP();
        }
        
        if(workflow.getInvertForVs() == 1){
            lastGradientVs = gradient.getVelocityS();
            lastConjugateGradientVs = gradient.getVelocityS();
        }
        
        if(workflow.getInvertForDensity() == 1){
            lastGradientDensity = gradient.getDensity();
            lastConjugateGradientDensity = gradient.getDensity();
        }
    
    }
    
    else {
    
        if(workflow.getInvertForVp() == 1){
            scai::lama::DenseVector<ValueType> gradientVp; 
            gradientVp = gradient.getVelocityP();
            scai::lama::DenseVector<ValueType> conjugateGradientVp;
            
            this->calcConjugateGradient(conjugateGradientVp, gradientVp, lastConjugateGradientVp, lastGradientVp);
            
            lastConjugateGradientVp = conjugateGradientVp;
            lastGradientVp = gradientVp;
            
            gradient.setVelocityP(conjugateGradientVp);
        }
        
        if(workflow.getInvertForVs() == 1){
            scai::lama::DenseVector<ValueType> gradientVs;
            gradientVs = gradient.getVelocityS();
            scai::lama::DenseVector<ValueType> conjugateGradientVs;
            
            this->calcConjugateGradient(conjugateGradientVs, gradientVs, lastConjugateGradientVs, lastGradientVs);
            
            lastConjugateGradientVs = conjugateGradientVs;
            lastGradientVs = gradientVs;
            
            gradient.setVelocityS(conjugateGradientVs);
        }
        
        if(workflow.getInvertForDensity() == 1){
            scai::lama::DenseVector<ValueType> gradientDensity;
            gradientDensity = gradient.getDensity();
            scai::lama::DenseVector<ValueType> conjugateGradientDensity;
            
            this->calcConjugateGradient(conjugateGradientDensity, gradientDensity, lastConjugateGradientDensity, lastGradientDensity);
            
            lastConjugateGradientVp = conjugateGradientDensity;
            lastGradientVp = gradientDensity;
            
            gradient.setDensity(conjugateGradientDensity);
        }
        
    }
    
    gradient.scale(model, workflow);
   
}


/*! \brief Calculate the conjugate gradient direction after Polak and Ribiere for one parameter
 * 
 * 
 \param conjugateGradient Dense Vector which serves as in- and output
 \param gradient Dense Vector
 \param lastConjugateGradient Dense Vector
 \param lastGradient Dense Vector
 */
template <typename ValueType>
void KITGPI::Optimization::ConjugateGradient<ValueType>::calcConjugateGradient(scai::lama::DenseVector<ValueType> &conjugateGradient, scai::lama::DenseVector<ValueType> const &gradient, scai::lama::DenseVector<ValueType> const &lastConjugateGradient, scai::lama::DenseVector<ValueType> const &lastGradient)
{      
    ValueType beta;
    
    beta = gradient.dotProduct(scai::lama::eval<scai::lama::DenseVector<ValueType>>(gradient - lastGradient));
    beta /= lastGradient.dotProduct(lastGradient);
    
    conjugateGradient = gradient + beta * lastConjugateGradient;
   
}

template class KITGPI::Optimization::ConjugateGradient<double>;
template class KITGPI::Optimization::ConjugateGradient<float>;
