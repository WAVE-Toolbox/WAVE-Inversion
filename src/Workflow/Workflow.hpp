#pragma once

#include <scai/lama.hpp>
#include <iostream>

#include <Common/HostPrint.hpp>
#include <Configuration/Configuration.hpp>
#include "../Misfit/Misfit.hpp"

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
            
            void init(KITGPI::Configuration::Configuration config);
            void changeStage(KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, ValueType &steplengthInit);
            void readFromFile(std::string workflowFilename);
            void printParameters(scai::dmemo::CommunicatorPtr comm) const;
            void printInvertForParameters(scai::dmemo::CommunicatorPtr comm) const;
            
            bool getInvertForVp() const;
            bool getInvertForVs() const;
            bool getInvertForDensity() const;
            bool getInvertForPorosity() const;
            bool getInvertForSaturation() const;
            bool getInvertForReflectivity() const;
            bool getInvertForSigma() const;
            bool getInvertForEpsilon() const;
            bool getInvertForTauSigma() const;
            bool getInvertForTauEpsilon() const;
            std::vector<bool> getInvertForParameters() const;
            void setInvertForParameters(std::vector<bool> setInvertForParameters);
            ValueType getRelativeMisfitChange() const;
            scai::IndexType getFilterOrder() const;
            ValueType getLowerCornerFreq() const;
            ValueType getUpperCornerFreq() const;
            ValueType getMinOffset() const;
            ValueType getMaxOffset() const;
            ValueType getTimeDampingFactor() const;
            scai::IndexType getSkipCount() const;
            scai::lama::DenseVector<ValueType> getWeightingFreq() const;
            scai::lama::DenseVector<ValueType> getFrequencyVec() const;
            
            KITGPI::Workflow::Workflow<ValueType> &operator=(KITGPI::Workflow::Workflow<ValueType> const &rhs);
            
            int maxStage;
            scai::IndexType workflowStage;
            scai::IndexType iteration;
            scai::IndexType skipCount;
            scai::IndexType skipDT;
            bool isSeismic;

        private:

            std::ifstream workflowFile;
            
            bool invertForVp = false;
            bool invertForVs = false;
            bool invertForDensity = false;
            bool invertForPorosity = false;
            bool invertForSaturation = false;
            bool invertForReflectivity = false;
            bool invertForSigma = false;
            bool invertForEpsilon = false;
            bool invertForTauSigma = false;
            bool invertForTauEpsilon = false;
            ValueType relativeMisfitChange; 
            scai::IndexType filterOrder;
            ValueType lowerCornerFreq;
            ValueType upperCornerFreq;
            ValueType minOffset;
            ValueType maxOffset;
            ValueType timeDampingFactor; 
            scai::lama::DenseVector<ValueType> weightingFreq; 
            scai::lama::DenseVector<ValueType> frequencyVector;           
        };
    }
}
