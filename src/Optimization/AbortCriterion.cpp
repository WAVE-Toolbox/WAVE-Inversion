#include "AbortCriterion.hpp"

/*! \brief Check abort criterion/criteria
 * 
 * 
 \param comm Communictaor 
 \param misfit reference of misfit is used because sometimes the misfitStorage has to be reset
 \param config configuration 
 \param steplengthInit reference of initial step lenght is used because it has ot be reset sometimes 
 \param workflow reference of workflow is used because sometimes private members have to be reset 
 \param workflowStage workflow stage number 
 \param iteration iteration number 
 */
template <typename ValueType>
bool KITGPI::AbortCriterion<ValueType>::check(scai::dmemo::CommunicatorPtr comm, KITGPI::Misfit::Misfit<ValueType> &misfit, KITGPI::Configuration::Configuration config, ValueType &steplengthInit, Workflow::Workflow<ValueType> &workflow, int workflowStage, int iteration)
{
    bool breakLoop = false;
    
    if ( (iteration > 1) && ( std::abs((misfit.getMisfitSum(iteration) - misfit.getMisfitSum(iteration - 2)))/(misfit.getMisfitSum(iteration - 2)) < workflow.relativeMisfitChange) ) 
    {        
        HOST_PRINT(comm, "\nAbort criterion fulfilled \n");
        HOST_PRINT(comm, "|Misfit(it)-Misfit(it-2)| / Misfit(it-2) < " << workflow.relativeMisfitChange << "\n");
        
        if(workflowStage != workflow.maxStage-1){
            HOST_PRINT(comm, "\nChange workflow stage\n");
            workflow.changeStage(config, misfit, steplengthInit);}
            
        breakLoop = true;
        
    }
    
    return breakLoop;
    
}


template class KITGPI::AbortCriterion<double>;
template class KITGPI::AbortCriterion<float>;
