#pragma once

#include <iostream>
#include <scai/lama.hpp>
#include <Common/HostPrint.hpp>
#include <Configuration/Configuration.hpp>

#include "Misfit.hpp"
#include "../Workflow/Workflow.hpp"
#include "../WorkflowEM/Workflow.hpp"

namespace KITGPI
{

    /*! \brief Class to check the abort criterion/criteria 
        *
        */
    
    template <typename ValueType>
    class AbortCriterion
    {

    public:

        /* Default constructor and destructor */
        AbortCriterion(){};
        ~AbortCriterion(){};
        
        bool check(scai::dmemo::CommunicatorPtr comm, KITGPI::Misfit::Misfit<ValueType> &misfit, KITGPI::Configuration::Configuration config, ValueType &steplengthInit, KITGPI::Workflow::Workflow<ValueType> &workflow, bool const breakLoopEM, scai::IndexType const breakLoopType);
        bool check(scai::dmemo::CommunicatorPtr comm, KITGPI::Misfit::Misfit<ValueType> &misfit, KITGPI::Configuration::Configuration configEM, ValueType &steplengthInitEM, Workflow::WorkflowEM<ValueType> &workflowEM, bool const breakLoop, scai::IndexType const breakLoopType);
        int const iterationStep = 2;
    };
}
