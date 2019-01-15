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
    if(!workflowFile){
        COMMON_THROWEXCEPTION("Workflow file not found! Abort program...")}
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
    workflowFile >> invertForVp >> invertForVs >> invertForDensity >> relativeMisfitChange >> filterOrder >> lowerCornerFreq >> upperCornerFreq >> useGradientTaper;
    
    workflowFile.close();
}

/*! \brief Print current set of workflow variables to terminal
 *
 \param comm Communicator 
 */
template <typename ValueType>
void KITGPI::Workflow::Workflow<ValueType>::printParameters(scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "=================================================\n");
    HOST_PRINT(comm, "Workflow stage parameters: \n");
    HOST_PRINT(comm, "invertForVp = " << invertForVp << "\n");
    HOST_PRINT(comm, "invertForVs = " << invertForVs << "\n");
    HOST_PRINT(comm, "invertForDensity = " << invertForDensity << "\n");
    HOST_PRINT(comm, "relativeMisfitChange = " << relativeMisfitChange << "\n");
    HOST_PRINT(comm, "filterOrder = " << filterOrder << "\n");
    HOST_PRINT(comm, "lowerCornerFreq = " << lowerCornerFreq << "\n");
    HOST_PRINT(comm, "upperCornerFreq = " << upperCornerFreq << "\n");
    HOST_PRINT(comm, "useGradientTaper = " << useGradientTaper << "\n");
}

/*! \brief Return copy of invertForVp
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForVp() const
{
    return invertForVp;
}

/*! \brief Return copy of invertForVs
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForVs() const
{
    return invertForVs;
}

/*! \brief Return copy of invertForDensity
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForDensity() const
{
    return invertForDensity;
}

/*! \brief Return copy of relativeMisfitChange
 */
template <typename ValueType>
ValueType KITGPI::Workflow::Workflow<ValueType>::getRelativeMisfitChange() const
{
    return relativeMisfitChange;
}

/*! \brief Return copy of filterOrder
 */
template <typename ValueType>
scai::IndexType KITGPI::Workflow::Workflow<ValueType>::getFilterOrder() const
{
    return filterOrder;
}

/*! \brief Return copy of lowerCornerFreq
 */
template <typename ValueType>
ValueType KITGPI::Workflow::Workflow<ValueType>::getLowerCornerFreq() const
{
    return lowerCornerFreq;
}

/*! \brief Return copy of upperCornerFreq
 */
template <typename ValueType>
ValueType KITGPI::Workflow::Workflow<ValueType>::getUpperCornerFreq() const
{
    return upperCornerFreq;
}

/*! \brief Return copy of invertForDensity
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getUseGradientTaper() const
{
    return useGradientTaper;
}
 
template class KITGPI::Workflow::Workflow<double>;
template class KITGPI::Workflow::Workflow<float>;
