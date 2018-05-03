#pragma once

#include <scai/lama.hpp>
#include <iostream>

#include <Configuration/Configuration.hpp>
#include <../Optimization/Misfit/Misfit.hpp>

namespace KITGPI
{
    
    //! \brief Workflow namespace
    namespace Workflow
    {
        /*! \brief Class to define and use workflow stages
         * 
         *  The workflow is read from a file (.txt) which has to be specified in the configuration file.
         *  It is possible, for example, to invert for different parameters in different workflow stages.
         */
        template <typename ValueType>
        class Workflow
        {

        public:
            
            /* Default constructor and destructor */
            Workflow(){};
            ~Workflow(){};

            Workflow(KITGPI::Configuration::Configuration config);
            
            void init(KITGPI::Configuration::Configuration config);
            void changeStage(KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, ValueType &steplengthInit);
            void readFromFile(std::string workflowFilename);
      
            int maxStage;
            
            bool invertForVp = false;
            bool invertForVs = false;
            bool invertForDensity = false;
            
        private:
      
            int currentStage;
            std::ifstream workflowFile;
   
            
        };
    }
}
