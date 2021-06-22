#pragma once

#include <scai/lama.hpp>
#include <iostream>

#include <Common/HostPrint.hpp>
#include <Configuration/Configuration.hpp>
#include <../Misfit/Misfit.hpp>

namespace KITGPI
{
    
    //! \brief Workflow namespace
    namespace Workflow
    {
        /*! \brief Class to define and use workflowEM stages
         * 
         *  The workflowEM is read from a file (.txt) which has to be specified in the configuration file.
         *  It is possible, for example, to invert for different parameters in different workflowEM stages.
         */
        template <typename ValueType>
        class WorkflowEM
        {

        public:
            
            /* Default constructor and destructor */
            WorkflowEM(){};
            ~WorkflowEM(){};

            WorkflowEM(KITGPI::Configuration::Configuration configEM);
            
            void init(KITGPI::Configuration::Configuration configEM);
            void changeStage(KITGPI::Configuration::Configuration configEM, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, ValueType &steplengthInitEM);
            void readFromFile(std::string workflowFilename);
            void printParameters(scai::dmemo::CommunicatorPtr comm);
            
            bool getInvertForSigmaEM() const;
            bool getInvertForEpsilonEM() const;
            bool getInvertForTauSigmaEM() const;
            bool getInvertForTauEpsilonEM() const;
            bool getInvertForPorosity() const;
            bool getInvertForSaturation() const;
            std::vector<bool> getInvertParameters() const;
            void setInvertParameters(std::vector<bool> setInvertParameters);
            ValueType getRelativeMisfitChange() const;
            scai::IndexType getFilterOrder() const;
            ValueType getLowerCornerFreq() const;
            ValueType getUpperCornerFreq() const;
            ValueType getTimeDampingFactor() const;
            scai::IndexType getSkipCount() const;
            
            int maxStage;
            scai::IndexType workflowStage;
            scai::IndexType iteration;
            scai::IndexType skipCount;

        private:

            std::ifstream workflowFile;
            
            bool invertForSigmaEM = true;
            bool invertForEpsilonEM = true;
            bool invertForTauSigmaEM = true;
            bool invertForTauEpsilonEM = true;
            bool invertForPorosity = true;
            bool invertForSaturation = true;
            ValueType relativeMisfitChange; 
            scai::IndexType filterOrder;
            ValueType lowerCornerFreq;
            ValueType upperCornerFreq;
            ValueType timeDampingFactor;
            
        };
    }
}
