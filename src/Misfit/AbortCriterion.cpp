#include "AbortCriterion.hpp"

/*! \brief Check abort criterion/criteria
 * 
 * 
 \param comm Communicator 
 \param misfit reference of misfit is used because sometimes the misfitStorage has to be reset
 \param config configuration 
 \param steplengthInit reference of initial step lenght is used because it has not be reset sometimes 
 \param workflow reference of workflow is used because sometimes private members have to be reset 
 */
template <typename ValueType>
bool KITGPI::AbortCriterion<ValueType>::check(scai::dmemo::CommunicatorPtr comm, KITGPI::Misfit::Misfit<ValueType> &misfit, KITGPI::Configuration::Configuration config, ValueType &steplengthInit, KITGPI::Workflow::Workflow<ValueType> &workflow, bool const breakLoopEM, scai::IndexType const breakLoopType)
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
        if (breakLoopType == 1 && breakLoopEM == true) {
            breakLoop = true;
        }
        if ((breakLoopType != 2 && breakLoop == true) || (breakLoopType == 2 && breakLoop == true && breakLoopEM == true)) {
            if(workflow.workflowStage != workflow.maxStage-1) {
                HOST_PRINT(comm, "\nChange workflow stage\n");
                workflow.changeStage(config, misfit, steplengthInit);                
            }                
        }        
    }    
    return breakLoop;    
}

/*! \brief Check abort criterion/criteria
 * 
 * 
 \param comm Communicator 
 \param misfit reference of misfit is used because sometimes the misfitStorage has to be reset
 \param configEM configuration 
 \param steplengthInitEM reference of initial step lenght is used because it has not be reset sometimes 
 \param workflowEM reference of workflowEM is used because sometimes private members have to be reset 
 */
template <typename ValueType>
bool KITGPI::AbortCriterion<ValueType>::check(scai::dmemo::CommunicatorPtr comm, KITGPI::Misfit::Misfit<ValueType> &misfit, KITGPI::Configuration::Configuration configEM, ValueType &steplengthInitEM, Workflow::WorkflowEM<ValueType> &workflowEM, bool const breakLoop, scai::IndexType const breakLoopType)
{
    bool breakLoopEM = false;
    
    if ( workflowEM.iteration > 1 ) {        
        if ( std::abs( misfit.getMisfitSum(workflowEM.iteration) - misfit.getMisfitSum(workflowEM.iteration-2) )/(misfit.getMisfitSum(workflowEM.iteration - 2)) < workflowEM.getRelativeMisfitChange()) {
            breakLoopEM = true;  
            HOST_PRINT(comm, "\nAbort criterion 1 fulfilled \n");
            HOST_PRINT(comm, "|Misfit(it)-Misfit(it-2)| / Misfit(it-2) < " << workflowEM.getRelativeMisfitChange() << "\n");
        }
        if ( workflowEM.iteration > iterationStep-1 && ( misfit.getMisfitSum(workflowEM.iteration) - misfit.getMisfitSum(workflowEM.iteration - iterationStep) )/(misfit.getMisfitSum(workflowEM.iteration - iterationStep)) > iterationStep * workflowEM.getRelativeMisfitChange()) {
            breakLoopEM = true;  
            HOST_PRINT(comm, "\nAbort criterion 2 fulfilled \n");
            HOST_PRINT(comm, "(Misfit(it)-Misfit(it-" << iterationStep << ")) / Misfit(it-"<< iterationStep << ") > " << iterationStep * workflowEM.getRelativeMisfitChange() << "\n");
        }        
        if (breakLoopType == 1 && breakLoop == true) {
            breakLoopEM = true;
        }
        if ((breakLoopType != 2 && breakLoopEM == true) || (breakLoopType == 2 && breakLoop == true && breakLoopEM == true)) {
            if(workflowEM.workflowStage != workflowEM.maxStage-1) {
                HOST_PRINT(comm, "\nChange workflowEM stage\n");
                workflowEM.changeStage(configEM, misfit, steplengthInitEM);                
            }                
        }      
    }    
    return breakLoopEM;    
}

template class KITGPI::AbortCriterion<double>;
template class KITGPI::AbortCriterion<float>;
