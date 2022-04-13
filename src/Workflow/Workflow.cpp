#include "Workflow.hpp"

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
    std::string equationType = config.get<std::string>("equationType");
    isSeismic = Common::checkEquationType<ValueType>(equationType);    
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
    skipCount = workflowStage;
    readFromFile(config.get<std::string>("workflowFilename"));
    
    dataMisfit.clearStorage();
    steplengthInit = config.get<ValueType>("steplengthInit");  
    
    IndexType gradientDomain = config.getAndCatch("gradientDomain", 0);
    IndexType dtinversion = config.get<IndexType>("DTInversion"); 
    if (dtinversion == 1 || upperCornerFreq == 0.0) {
        skipDT = 1;
    } else {
        ValueType dtMax = 1.0 / (upperCornerFreq * 5);
        skipDT = floor(dtMax / config.get<ValueType>("DT")); // Nyquist sampling step
        if (gradientDomain != 0) {
            scai::IndexType nFFT = Common::calcNextPowTwo<ValueType>(skipDT - 1);
            if (nFFT > skipDT) {
                skipDT = nFFT / 2;
            } else {
                skipDT = nFFT;
            }
        }
    }
        
    IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    if (gradientDomain != 0) {
        IndexType NT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
        scai::IndexType nFFT = Common::calcNextPowTwo<ValueType>(NT - 1);
        SCAI_ASSERT_ERROR(nFFT == NT, "NT must be the power of 2 when gradientDomain != 0!");
        ValueType df = 1 / (nFFT * config.get<ValueType>("DT"));
        IndexType fc1Ind = ceil(lowerCornerFreq / df);
        IndexType fc2Ind = ceil(upperCornerFreq / df);  
        if (useSourceEncode == 0) {
            fc1Ind = ceil((ValueType)fc1Ind / 2); // to cover all effective frequencies
            fc2Ind *= 2; // to cover all effective frequencies
            if (fc1Ind == 0)
                fc1Ind = 1;
            IndexType nfc12 = ceil(ValueType(fc2Ind-fc1Ind+1));  
            frequencyVector = scai::lama::linearDenseVector<ValueType>(nfc12, fc1Ind*df, df);
            scai::lama::DenseVector<ValueType> temp = scai::lama::linearDenseVector<ValueType>(fc1Ind, 0, df);
            frequencyVector.cat(temp, frequencyVector);
            weightingFreq = lama::linearDenseVector<ValueType>(nfc12, 1, 0);
            temp = scai::lama::linearDenseVector<ValueType>(fc1Ind, 0, 0);
            weightingFreq.cat(temp, weightingFreq);
        } else {
            if (fc1Ind == 0)
                fc1Ind = 2; // ensure steady-state wavefields.
            if (fc1Ind % 2 == 1)
                fc1Ind += 1; // ensure steady-state wavefields.
            IndexType nfc12 = ceil(ValueType(fc2Ind-fc1Ind+1)/2);  
            frequencyVector = scai::lama::linearDenseVector<ValueType>(nfc12, fc1Ind*df, 2*df);
            scai::lama::DenseVector<ValueType> temp = scai::lama::linearDenseVector<ValueType>(fc1Ind/2, 0, 2*df);
            frequencyVector.cat(temp, frequencyVector);
            ValueType average = (1.0 + 1.0 / maxStage) / 2;
            ValueType start = 1.0 - ValueType(workflowStage) / maxStage;    
            ValueType step = (average - start) / (nfc12 - 1);
            weightingFreq = lama::linearDenseVector<ValueType>(nfc12, start, 2*step);
            weightingFreq = 1;
            temp = scai::lama::linearDenseVector<ValueType>(fc1Ind/2, 0, 0);
            weightingFreq.cat(temp, weightingFreq);
        }
    }
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
        workflowFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');        
    }
        
    /* Extract variables from current workflow stage */
    if (isSeismic)
        workflowFile >> invertForVp >> invertForVs >> invertForDensity >> invertForPorosity >> invertForSaturation >> relativeMisfitChange >> filterOrder >> lowerCornerFreq >> upperCornerFreq >> timeDampingFactor;
    else
        workflowFile >> invertForSigmaEM >> invertForEpsilonEM >> invertForTauSigmaEM >> invertForTauEpsilonEM >> invertForPorosity >> invertForSaturation >> relativeMisfitChange >> filterOrder >> lowerCornerFreq >> upperCornerFreq >> timeDampingFactor;
    
    workflowFile.close();
}

/*! \brief Print current set of workflow variables to terminal
 *
 \param comm Communicator 
 */
template <typename ValueType>
void KITGPI::Workflow::Workflow<ValueType>::printInvertForParameters(scai::dmemo::CommunicatorPtr comm) const
{
    if (isSeismic) {
        if (invertForVp)
            HOST_PRINT(comm, "invertForVp = " << invertForVp << "\n");
        if (invertForVs)
            HOST_PRINT(comm, "invertForVs = " << invertForVs << "\n");
        if (invertForDensity)
            HOST_PRINT(comm, "invertForDensity = " << invertForDensity << "\n");
    } else {
        if (invertForSigmaEM)
            HOST_PRINT(comm, "invertForSigmaEM = " << invertForSigmaEM << "\n");
        if (invertForEpsilonEM)
            HOST_PRINT(comm, "invertForEpsilonEM = " << invertForEpsilonEM << "\n");
        if (invertForTauSigmaEM)
            HOST_PRINT(comm, "invertForTauSigmaEM = " << invertForTauSigmaEM << "\n");
        if (invertForTauEpsilonEM)
            HOST_PRINT(comm, "invertForTauEpsilonEM = " << invertForTauEpsilonEM << "\n");
    }
    if (invertForPorosity)
        HOST_PRINT(comm, "invertForPorosity = " << invertForPorosity << "\n");
    if (invertForSaturation)
        HOST_PRINT(comm, "invertForSaturation = " << invertForSaturation << "\n");
    if (invertForReflectivity)
        HOST_PRINT(comm, "invertForReflectivity = " << invertForReflectivity << "\n");
}

/*! \brief Print current set of workflow variables to terminal
 *
 \param comm Communicator 
 */
template <typename ValueType>
void KITGPI::Workflow::Workflow<ValueType>::printParameters(scai::dmemo::CommunicatorPtr comm) const
{
    HOST_PRINT(comm, "=================================================\n");
    HOST_PRINT(comm, "Workflow stage parameters: \n");
    printInvertForParameters(comm);
    HOST_PRINT(comm, "relativeMisfitChange = " << relativeMisfitChange << "\n");
    HOST_PRINT(comm, "filterOrder = " << filterOrder << "\n");
    HOST_PRINT(comm, "lowerCornerFreq = " << lowerCornerFreq << "\n");
    HOST_PRINT(comm, "upperCornerFreq = " << upperCornerFreq << "\n");
    HOST_PRINT(comm, "timeDampingFactor = " << timeDampingFactor << "\n");
    if (weightingFreq.size() != 0) {
        HOST_PRINT(comm, "frequencyVector =");
        for (int i=0; i<frequencyVector.size(); i++) {
            HOST_PRINT(comm, " " << frequencyVector[i]);
        }
        HOST_PRINT(comm, "\n");
        HOST_PRINT(comm, "weightingFreq =");
        for (int i=0; i<weightingFreq.size(); i++) {
            HOST_PRINT(comm, " " << weightingFreq[i]);
        }
        HOST_PRINT(comm, "\n");
    }
}

/*! \brief Return copy of invertForVp
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForVp() const
{
    SCAI_ASSERT_ERROR(isSeismic, "!isSeismic");
    return invertForVp;
}

/*! \brief Return copy of invertForVs
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForVs() const
{
    SCAI_ASSERT_ERROR(isSeismic, "!isSeismic");
    return invertForVs;
}

/*! \brief Return copy of invertForDensity
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForDensity() const
{
    SCAI_ASSERT_ERROR(isSeismic, "!isSeismic");
    return invertForDensity;
}

/*! \brief Return copy of invertForPorosity
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForPorosity() const
{
    return invertForPorosity;
}

/*! \brief Return copy of invertForSaturation
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForSaturation() const
{
    return invertForSaturation;
}

/*! \brief Return copy of invertForReflectivity
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForReflectivity() const
{
    return invertForReflectivity;
}

/*! \brief Return copy of invertForSigmaEM
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForSigmaEM() const
{
    SCAI_ASSERT_ERROR(!isSeismic, "isSeismic");
    return invertForSigmaEM;
}

/*! \brief Return copy of invertForEpsilonEM
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForEpsilonEM() const
{
    SCAI_ASSERT_ERROR(!isSeismic, "isSeismic");
    return invertForEpsilonEM;
}

/*! \brief Return copy of invertForTauSigmaEM
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForTauSigmaEM() const
{
    SCAI_ASSERT_ERROR(!isSeismic, "isSeismic");
    return invertForTauSigmaEM;
}

/*! \brief Return copy of invertForTauEpsilonEM
 */
template <typename ValueType>
bool KITGPI::Workflow::Workflow<ValueType>::getInvertForTauEpsilonEM() const
{
    return invertForTauEpsilonEM;
}

/*! \brief Return the vector of the inverted parameters
 */
template <typename ValueType>
std::vector<bool> KITGPI::Workflow::Workflow<ValueType>::getInvertForParameters() const
{
    std::vector<bool> invertForParameters;
    if (isSeismic) {
        std::vector<bool> temp{invertForVp, invertForVs, invertForDensity, invertForPorosity, invertForSaturation, invertForReflectivity};
        invertForParameters = temp;
    } else {
        std::vector<bool> temp{invertForSigmaEM, invertForEpsilonEM, invertForTauSigmaEM, invertForTauEpsilonEM, invertForPorosity, invertForSaturation, invertForReflectivity};
        invertForParameters = temp;
    }
    return invertForParameters;
}

/*! \brief Set the vector of the inverted parameters
 */
template <typename ValueType>
void KITGPI::Workflow::Workflow<ValueType>::setInvertForParameters(std::vector<bool> setInvertForParameters)
{
    if (isSeismic) {
        invertForVp = setInvertForParameters[0];
        invertForVs = setInvertForParameters[1];
        invertForDensity = setInvertForParameters[2];
        invertForPorosity = setInvertForParameters[3];
        invertForSaturation = setInvertForParameters[4];
        invertForReflectivity = setInvertForParameters[5];
    } else {    
        invertForSigmaEM = setInvertForParameters[0];
        invertForEpsilonEM = setInvertForParameters[1];
        invertForTauSigmaEM = setInvertForParameters[2];
        invertForTauEpsilonEM = setInvertForParameters[3];
        invertForPorosity = setInvertForParameters[4];
        invertForSaturation = setInvertForParameters[5];
        invertForReflectivity = setInvertForParameters[6];
    }
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

/*! \brief Return copy of timeDampingFactor
 */
template <typename ValueType>
ValueType KITGPI::Workflow::Workflow<ValueType>::getTimeDampingFactor() const
{
    return timeDampingFactor;
}

/*! \brief Return copy of skipCount
 */
template <typename ValueType>
scai::IndexType KITGPI::Workflow::Workflow<ValueType>::getSkipCount() const
{
    return skipCount;
}

/*! \brief Return copy of weightingFreq
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Workflow::Workflow<ValueType>::getWeightingFreq() const
{
    return weightingFreq;
}

/*! \brief Return copy of frequencyVector
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Workflow::Workflow<ValueType>::getFrequencyVec() const
{
    return frequencyVector;
}

/*! \brief Overloading = Operation
 *
 \param rhs Workflow which is copied.
 */
template <typename ValueType>
KITGPI::Workflow::Workflow<ValueType> &KITGPI::Workflow::Workflow<ValueType>::operator=(KITGPI::Workflow::Workflow<ValueType> const &rhs)
{
    invertForVp = rhs.invertForVp;
    invertForVs = rhs.invertForVs;
    invertForDensity = rhs.invertForDensity;
    invertForPorosity = rhs.invertForPorosity;
    invertForSaturation = rhs.invertForSaturation;
    invertForReflectivity = rhs.invertForReflectivity;
    invertForSigmaEM = rhs.invertForSigmaEM;
    invertForEpsilonEM = rhs.invertForEpsilonEM;
    invertForTauSigmaEM = rhs.invertForTauSigmaEM;
    invertForTauEpsilonEM = rhs.invertForTauEpsilonEM;
    
    return *this;
}

template class KITGPI::Workflow::Workflow<double>;
template class KITGPI::Workflow::Workflow<float>;
