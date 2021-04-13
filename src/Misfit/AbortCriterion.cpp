#include "AbortCriterion.hpp"

/*! \brief Check abort criterion/criteria
 * 
 * 
 \param comm Communicator 
 \param misfit reference of misfit is used because sometimes the misfitStorage has to be reset
 \param config configuration 
 \param steplengthInit reference of initial step lenght is used because it has ot be reset sometimes 
 \param workflow reference of workflow is used because sometimes private members have to be reset 
 */
template <typename ValueType>
bool KITGPI::AbortCriterion<ValueType>::check(scai::dmemo::CommunicatorPtr comm, KITGPI::Misfit::Misfit<ValueType> &misfit, KITGPI::Configuration::Configuration config, ValueType &steplengthInit, Workflow::Workflow<ValueType> &workflow)
{
    bool breakLoop = false;
    
    if ( workflow.iteration > 1 ) {        
        if ( std::abs( misfit.getMisfitSum(workflow.iteration) - misfit.getMisfitSum(workflow.iteration-2) )/(misfit.getMisfitSum(workflow.iteration - 2)) < workflow.getRelativeMisfitChange()) {
            breakLoop = true;  
            HOST_PRINT(comm, "\nAbort criterion 1 fulfilled \n");
            HOST_PRINT(comm, "|Misfit(it)-Misfit(it-2)| / Misfit(it-2) < " << workflow.getRelativeMisfitChange() << "\n");
        }
        if ( workflow.iteration > iterationStep-1 && ( misfit.getMisfitSum(workflow.iteration) - misfit.getMisfitSum(workflow.iteration - iterationStep) )/(misfit.getMisfitSum(workflow.iteration - iterationStep)) > iterationStep * workflow.getRelativeMisfitChange()) {
            breakLoop = true;  
            HOST_PRINT(comm, "\nAbort criterion 2 fulfilled \n");
            HOST_PRINT(comm, "(Misfit(it)-Misfit(it-" << iterationStep << ")) / Misfit(it-"<< iterationStep << ") > " << iterationStep * workflow.getRelativeMisfitChange() << "\n");
        }    
        if(workflow.workflowStage != workflow.maxStage-1){
            HOST_PRINT(comm, "\nChange workflow stage\n");
            workflow.changeStage(config, misfit, steplengthInit);}
            
        breakLoop = true;        
    }    
    return breakLoop;    
}

template class KITGPI::AbortCriterion<double>;
template class KITGPI::AbortCriterion<float>;
