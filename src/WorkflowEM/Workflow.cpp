#include "Workflow.hpp"

/*! \brief Constructor
 *
 *
 \param configEM Configuration
 */
template <typename ValueType>
KITGPI::Workflow::WorkflowEM<ValueType>::WorkflowEM(KITGPI::Configuration::Configuration configEM)
{
    init(configEM);
}

/*! \brief Initialize members of the object
 *
 *
 \param configEM Configuration
 */
template <typename ValueType>
void KITGPI::Workflow::WorkflowEM<ValueType>::init(KITGPI::Configuration::Configuration configEM)
{
    workflowFile.open(configEM.get<std::string>("workflowFilename"));
    if(!workflowFile){
        COMMON_THROWEXCEPTION("Workflow file not found! Abort program...")}
    maxStage = std::count(std::istreambuf_iterator<char>(workflowFile), std::istreambuf_iterator<char>(), '\n');
    maxStage = maxStage-2; // substract comment lines
    workflowFile.close();
    skipCount = 0;
    readFromFile(configEM.get<std::string>("workflowFilename"));
    
}

/*! \brief Reset and adapt parameters when workflowEM is changed
 *
 * 
 \param configEM Configuration
 \param dataMisfitEM data misfit
 \param steplengthInitEM initial step length specified in the configuration file
 */
template <typename ValueType>
void KITGPI::Workflow::WorkflowEM<ValueType>::changeStage(KITGPI::Configuration::Configuration configEM, KITGPI::Misfit::Misfit<ValueType> &dataMisfitEM, ValueType &steplengthInitEM)
{
    skipCount += 1;
    readFromFile(configEM.get<std::string>("workflowFilename"));
    
    dataMisfitEM.clearStorage();
    steplengthInitEM = configEM.get<ValueType>("steplengthInit");  
}

/*! \brief Read parameters from workflowEM file
 *
 * 
 \param workflowFilename name of the workflowEM file specified in the configuration file
 */
template <typename ValueType>
void KITGPI::Workflow::WorkflowEM<ValueType>::readFromFile(std::string workflowFilename)
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
        
    /* Extract variables from current workflowEM stage */
    workflowFile >> invertForSigmaEM >> invertForEpsilonEM >> invertForTauSigmaEM >> invertForTauEpsilonEM >> invertForPorosity >> invertForSaturation >> relativeMisfitChange >> filterOrder >> lowerCornerFreq >> upperCornerFreq >> timeDampingFactor;
    
    workflowFile.close();
}

/*! \brief Print current set of workflowEM variables to terminal
 *
 \param comm Communicator 
 */
template <typename ValueType>
void KITGPI::Workflow::WorkflowEM<ValueType>::printParameters(scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "=================================================\n");
    HOST_PRINT(comm, "Workflow stage parameters: \n");
    HOST_PRINT(comm, "invertForSigmaEM = " << invertForSigmaEM << "\n");
    HOST_PRINT(comm, "invertForEpsilonEM = " << invertForEpsilonEM << "\n");
    HOST_PRINT(comm, "invertForTauSigmaEM = " << invertForTauSigmaEM << "\n");
    HOST_PRINT(comm, "invertForTauEpsilonEM = " << invertForTauEpsilonEM << "\n");
    HOST_PRINT(comm, "invertForPorosity = " << invertForPorosity << "\n");
    HOST_PRINT(comm, "invertForSaturation = " << invertForSaturation << "\n");
    HOST_PRINT(comm, "relativeMisfitChange = " << relativeMisfitChange << "\n");
    HOST_PRINT(comm, "filterOrder = " << filterOrder << "\n");
    HOST_PRINT(comm, "lowerCornerFreq = " << lowerCornerFreq << "\n");
    HOST_PRINT(comm, "upperCornerFreq = " << upperCornerFreq << "\n");
    HOST_PRINT(comm, "timeDampingFactor = " << timeDampingFactor << "\n");
}

/*! \brief Return copy of invertForSigmaEM
 */
template <typename ValueType>
bool KITGPI::Workflow::WorkflowEM<ValueType>::getInvertForSigmaEM() const
{
    return invertForSigmaEM;
}

/*! \brief Return copy of invertForEpsilonEM
 */
template <typename ValueType>
bool KITGPI::Workflow::WorkflowEM<ValueType>::getInvertForEpsilonEM() const
{
    return invertForEpsilonEM;
}

/*! \brief Return copy of invertForTauSigmaEM
 */
template <typename ValueType>
bool KITGPI::Workflow::WorkflowEM<ValueType>::getInvertForTauSigmaEM() const
{
    return invertForTauSigmaEM;
}

/*! \brief Return copy of invertForTauEpsilonEM
 */
template <typename ValueType>
bool KITGPI::Workflow::WorkflowEM<ValueType>::getInvertForTauEpsilonEM() const
{
    return invertForTauEpsilonEM;
}

/*! \brief Return copy of invertForPorosity
 */
template <typename ValueType>
bool KITGPI::Workflow::WorkflowEM<ValueType>::getInvertForPorosity() const
{
    return invertForPorosity;
}

/*! \brief Return copy of invertForSaturation
 */
template <typename ValueType>
bool KITGPI::Workflow::WorkflowEM<ValueType>::getInvertForSaturation() const
{
    return invertForSaturation;
}

/*! \brief Return the vector of the inverted parameters
 */
template <typename ValueType>
std::vector<bool> KITGPI::Workflow::WorkflowEM<ValueType>::getInvertParameters() const
{
    std::vector<bool> invertParameters{invertForSigmaEM, invertForEpsilonEM, invertForTauSigmaEM, invertForTauEpsilonEM, invertForPorosity, invertForSaturation};
    return invertParameters;
}

/*! \brief Set the vector of the inverted parameters
 */
template <typename ValueType>
void KITGPI::Workflow::WorkflowEM<ValueType>::setInvertParameters(std::vector<bool> setInvertParameters)
{
    invertForSigmaEM = setInvertParameters[0];
    invertForEpsilonEM = setInvertParameters[1];
    invertForTauSigmaEM = setInvertParameters[2];
    invertForTauEpsilonEM = setInvertParameters[3];
    invertForPorosity = setInvertParameters[4];
    invertForSaturation = setInvertParameters[5];
}

/*! \brief Return copy of relativeMisfitChange
 */
template <typename ValueType>
ValueType KITGPI::Workflow::WorkflowEM<ValueType>::getRelativeMisfitChange() const
{
    return relativeMisfitChange;
}

/*! \brief Return copy of filterOrder
 */
template <typename ValueType>
scai::IndexType KITGPI::Workflow::WorkflowEM<ValueType>::getFilterOrder() const
{
    return filterOrder;
}

/*! \brief Return copy of lowerCornerFreq
 */
template <typename ValueType>
ValueType KITGPI::Workflow::WorkflowEM<ValueType>::getLowerCornerFreq() const
{
    return lowerCornerFreq;
}

/*! \brief Return copy of upperCornerFreq
 */
template <typename ValueType>
ValueType KITGPI::Workflow::WorkflowEM<ValueType>::getUpperCornerFreq() const
{
    return upperCornerFreq;
}

/*! \brief Return copy of timeDampingFactor
 */
template <typename ValueType>
ValueType KITGPI::Workflow::WorkflowEM<ValueType>::getTimeDampingFactor() const
{
    return timeDampingFactor;
}

/*! \brief Return copy of skipCount
 */
template <typename ValueType>
scai::IndexType KITGPI::Workflow::WorkflowEM<ValueType>::getSkipCount() const
{
    return skipCount;
}

template class KITGPI::Workflow::WorkflowEM<double>;
template class KITGPI::Workflow::WorkflowEM<float>;
