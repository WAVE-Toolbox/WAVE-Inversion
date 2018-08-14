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
            void printParameters(scai::dmemo::CommunicatorPtr comm);
            
            bool getInvertForVp() const;
            bool getInvertForVs() const;
            bool getInvertForDensity() const;
            ValueType getRelativeMisfitChange() const;
            scai::IndexType getFilterOrder() const;
            ValueType getLowerCornerFreq() const;
            ValueType getUpperCornerFreq() const;
            
            int maxStage;
            scai::IndexType workflowStage;
            scai::IndexType iteration;
            scai::IndexType skipCount;

        private:

            std::ifstream workflowFile;
            
            bool invertForVp;
            bool invertForVs;
            bool invertForDensity;
            ValueType relativeMisfitChange; 
            scai::IndexType filterOrder;
            ValueType lowerCornerFreq;
            ValueType upperCornerFreq;
            
        };
    }
}
