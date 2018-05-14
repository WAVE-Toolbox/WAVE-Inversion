#pragma once

#include <iostream>
#include <scai/lama.hpp>
#include <Common/HostPrint.hpp>
#include <Configuration/Configuration.hpp>

#include "Misfit/Misfit.hpp"
#include "../Workflow/Workflow.hpp"

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
        
        bool check(scai::dmemo::CommunicatorPtr comm, KITGPI::Misfit::Misfit<ValueType> &misfit, KITGPI::Configuration::Configuration config, ValueType &steplengthInit, Workflow::Workflow<ValueType> &workflow, int workflowStage, int iteration);


    };
}
