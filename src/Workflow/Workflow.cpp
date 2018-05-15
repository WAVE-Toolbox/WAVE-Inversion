#include "Workflow.hpp"

/*! \brief Constructor
 *
 *
 \param config Configuration
 */
template <typename ValueType>
KITGPI::Workflow::Workflow<ValueType>::Workflow(KITGPI::Configuration::Configuration config)
{
    init(config);
}

/*! \brief Initialize members of the object
 *
 *
 \param config Configuration
 */
template <typename ValueType>
void KITGPI::Workflow::Workflow<ValueType>::init(KITGPI::Configuration::Configuration config)
{
    workflowFile.open(config.get<std::string>("workflowFilename"));
    maxStage = std::count(std::istreambuf_iterator<char>(workflowFile), std::istreambuf_iterator<char>(), '\n');
    maxStage = maxStage-2; // substract comment lines
    workflowFile.close();
    skipCount = 0;
    readFromFile(config.get<std::string>("workflowFilename"));
    
}

/*! \brief Reset and adapt parameters when workflow is changed
 *
 * 
 \param config Configuration
 \param dataMisfit data misfit
 \param steplengthInit initial step length specified in the configuration file
 */
template <typename ValueType>
void KITGPI::Workflow::Workflow<ValueType>::changeStage(KITGPI::Configuration::Configuration config, KITGPI::Misfit::Misfit<ValueType> &dataMisfit, ValueType &steplengthInit)
{
    skipCount += 1;
    readFromFile(config.get<std::string>("workflowFilename"));
    
    dataMisfit.clearStorage();
    steplengthInit = config.get<ValueType>("steplengthInit");  
}

/*! \brief Read parameters from workflow file
 *
 * 
 \param workflowFilename name of the workflow file specified in the configuration file
 */
template <typename ValueType>
void KITGPI::Workflow::Workflow<ValueType>::readFromFile(std::string workflowFilename)
{   
    workflowFile.open(workflowFilename);
    
    /* Skip comment lines */
    char firstChar = workflowFile.peek();
    while (firstChar == '#') {
        workflowFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        firstChar = workflowFile.peek();
    }
    
    /* Skip stage lines (except for first stage) */
    for (int i = 0; i < skipCount; i++){
        workflowFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');}
        
    /* Extract variables from current workflow stage */
    workflowFile >> invertForVp >> invertForVs >> invertForDensity >> relativeMisfitChange;
    
    workflowFile.close();
}


template class KITGPI::Workflow::Workflow<double>;
template class KITGPI::Workflow::Workflow<float>;
