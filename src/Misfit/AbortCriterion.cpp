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
bool KITGPI::AbortCriterion<ValueType>::check(scai::dmemo::CommunicatorPtr comm, KITGPI::Misfit::Misfit<ValueType> &misfit, ValueType &steplengthInit, KITGPI::Workflow::Workflow<ValueType> &workflow, bool const breakLoopEM, scai::IndexType const breakLoopType)
{
    bool breakLoop = false;   
    if (workflow.iteration > 1) { 
        ValueType misfitResidual1 = 0;   
        ValueType misfitResidual2 = 0;      
        
        if (misfit.saveMultiMisfits || misfit.misfitType.length() > 2) {
            misfitResidual1 = misfit.getMisfitResidualMax(workflow.iteration - 1, workflow.iteration);
            misfitResidual2 = misfit.getMisfitResidualMax(workflow.iteration - 2, workflow.iteration);
        } else {
            misfitResidual1 = (misfit.getMisfitSum(workflow.iteration - 1) - misfit.getMisfitSum(workflow.iteration))/misfit.getMisfitSum(workflow.iteration - 1);
            misfitResidual2 = (misfit.getMisfitSum(workflow.iteration - 2) - misfit.getMisfitSum(workflow.iteration))/misfit.getMisfitSum(workflow.iteration - 2);
        }
    
        if ( std::abs(misfitResidual1) < workflow.getRelativeMisfitChange() && std::abs(misfitResidual2) < workflow.getRelativeMisfitChange()) {
            breakLoop = true;  
            HOST_PRINT(comm, "\nAbort criterion 1 fulfilled \n");
            HOST_PRINT(comm, "|Misfit(it-2)-Misfit(it)| / Misfit(it-2) < " << workflow.getRelativeMisfitChange() << "\n");
        }
        misfitResidual1 = -misfitResidual1;
        misfitResidual2 = -misfitResidual2;
        if (misfitResidual1 > 1.5 * workflow.getRelativeMisfitChange() && misfitResidual2 > 1.5 * workflow.getRelativeMisfitChange()) {
            breakLoop = true;  
            HOST_PRINT(comm, "\nAbort criterion 2 fulfilled \n");
            HOST_PRINT(comm, "(Misfit(it)-Misfit(it-2)) / Misfit(it-2) > " << 1.5 * workflow.getRelativeMisfitChange() << "\n");
        }        
        if (breakLoopType == 1 && breakLoopEM == true) {
            breakLoop = true;
        }      
    }    
    return breakLoop;    
}

template class KITGPI::AbortCriterion<double>;
template class KITGPI::AbortCriterion<float>;
