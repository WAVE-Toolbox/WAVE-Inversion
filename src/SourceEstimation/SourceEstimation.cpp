#include "SourceEstimation.hpp"

using namespace scai;

/*! \brief Initialize based on config which also initializes sourceSignalTaper
 \param config Configuration
 \param ctx Context Ptr
 \param sourceDistribution Source distribution pointer
 \param sourceSignalTaper 1D TaperPtr
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::init(Configuration::Configuration const &config, KITGPI::Workflow::Workflow<ValueType> const &workflow, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr sourceDistribution, Taper::Taper1D<ValueType> &sourceSignalTaper)
{
    std::string equationType = config.get<std::string>("equationType"); 
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);  
    isSeismic = Common::checkEquationType<ValueType>(equationType);
        
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    if (config.get<IndexType>("useSeismogramTaper") != 0 && config.get<IndexType>("useSeismogramTaper") != 2 && config.get<IndexType>("useSeismogramTaper") != 5) {
        if (config.get<IndexType>("useSeismogramTaper") == 4) {        
            init(tStepEnd, sourceDistribution, config.get<ValueType>("waterLevel"), config.get<std::string>("seismogramTaperName") + ".SrcEst");
        } else {     
            init(tStepEnd, sourceDistribution, config.get<ValueType>("waterLevel"), config.get<std::string>("seismogramTaperName"));
        }
    } else {
        init(tStepEnd, sourceDistribution, config.get<ValueType>("waterLevel"));
    }
    
    if (workflow.getMinOffset() != 0 || workflow.getMaxOffset() != 0)
        useOffsetMutes = true;
    
    useSourceEncode = config.getAndCatch("useSourceEncode", 0);

    if (config.get<IndexType>("useSourceSignalTaper") == 1) {
        sourceSignalTaper.init(std::make_shared<dmemo::NoDistribution>(tStepEnd), ctx, 1);
        if (config.get<IndexType>("sourceSignalTaperStart2") == 0 && config.get<IndexType>("sourceSignalTaperEnd2") == 0)
            sourceSignalTaper.calcCosineTaper(config.get<IndexType>("sourceSignalTaperStart1"), config.get<IndexType>("sourceSignalTaperEnd1"), 0);
        else
            sourceSignalTaper.calcCosineTaper(config.get<IndexType>("sourceSignalTaperStart1"), config.get<IndexType>("sourceSignalTaperEnd1"), config.get<IndexType>("sourceSignalTaperStart2"), config.get<IndexType>("sourceSignalTaperEnd2"), 0);
    }
}

/*! \brief Initialize based on raw variables
 \param nt Number of time steps
 \param sourceDistribution Source distribution pointer
 \param waterLvl Water level
 \param tprName if true a taper 
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::init(IndexType nt, dmemo::DistributionPtr sourceDistribution, ValueType waterLvl, std::string tprName)
{
    nFFT = Common::calcNextPowTwo<ValueType>(nt - 1);
    waterLevel = waterLvl;
    filter.allocate(sourceDistribution, std::make_shared<dmemo::NoDistribution>(nFFT));
    IndexType ns = sourceDistribution->getGlobalSize();
    refTraces.clear();
    mutes.clear();
    offsets.clear();
    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        scai::lama::DenseMatrix<ValueType> refTracesTemp(sourceDistribution, std::make_shared<dmemo::NoDistribution>(nt));
        for (scai::IndexType shotInd = 0; shotInd < ns; shotInd++) {
            scai::lama::DenseVector<ValueType> temp(nt, 0);
            refTracesTemp.setRow(temp, shotInd, common::BinaryOp::COPY);
        }
        refTraces.push_back(refTracesTemp);
    }
    for (scai::IndexType shotInd = 0; shotInd < ns; shotInd++) {
        scai::lama::DenseVector<ValueType> temp;
        mutes.push_back(temp);
        offsets.push_back(temp);
    }
    if (!tprName.empty()) {
        readTaper = true;
        taperName = tprName;
    }
}

/*! \brief Calculate the Wiener filter
 \param receivers Synthetic receivers
 \param receiversTrue Observed receivers
 \param shotNumber Shot number of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::estimateSourceSignal(KITGPI::Acquisition::Receivers<ValueType> const &receivers, KITGPI::Acquisition::Receivers<ValueType> const &receiversTrue, scai::IndexType shotInd, scai::IndexType shotNumber)
{
    lama::DenseVector<ComplexValueType> filterTmp1(filter.getColDistributionPtr(), 0.0);
    lama::DenseVector<ComplexValueType> filterTmp2(filter.getColDistributionPtr(), 0.0);

    addComponents(filterTmp1, receivers, receivers, shotInd, shotNumber);
//     filterTmp1 += waterLevel * receivers.getSeismogramHandler().getNumTracesTotal();
    filterTmp1 += waterLevel * filterTmp1.maxNorm();
    filterTmp1.unaryOp(filterTmp1, common::UnaryOp::RECIPROCAL);
    addComponents(filterTmp2, receivers, receiversTrue, shotInd, shotNumber);
    filterTmp1 *= filterTmp2;
    
    filter.setRow(filterTmp1, shotInd, common::BinaryOp::COPY);
}

/*! \brief Calculate the Wiener filter
 \param receivers Synthetic receivers
 \param receiversTrue Observed receivers
 \param shotNumber Shot number of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::estimateSourceSignalEncode(scai::dmemo::CommunicatorPtr commShot, scai::IndexType shotNumberEncode, KITGPI::Configuration::Configuration const &config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, KITGPI::Acquisition::Receivers<ValueType> const &receivers, KITGPI::Acquisition::Receivers<ValueType> const &receiversTrue)
{
    SCAI_ASSERT_ERROR(config.get<scai::IndexType>("useSourceSignalTaper") == 0, "useSourceSignalTaper != 0"); 
    if (useSourceEncode != 0) {
        double start_t_shot, end_t_shot; /* For timing */
        start_t_shot = common::Walltime::get();
        
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings;  
        Acquisition::Sources<ValueType> sources;
        ValueType shotIncr = 0; 
        sources.getAcquisitionSettings(config, shotIncr); // to get numshots
        sourceSettings = sources.getSourceSettings();
        std::vector<scai::IndexType> uniqueShotNos;
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
        IndexType numshotsIncr = sourceSettingsEncode.size();
        KITGPI::Acquisition::Receivers<ValueType> receiversSyn;
        KITGPI::Acquisition::Receivers<ValueType> receiversObs;
        IndexType countDecode = 0;
        for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
            if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumberEncode) {                
                for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
                    if (receivers.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
                        receiversSyn.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData() = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getDataDecode(countDecode);
                        receiversObs.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData() = receiversTrue.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getDataDecode(countDecode);
                    }
                }
                IndexType shotNumber = uniqueShotNos[sourceSettingsEncode[shotInd].row]; 
                                              
                estimateSourceSignal(receiversSyn, receiversObs, shotInd, shotNumber);
                countDecode++;
            }
        }
        end_t_shot = common::Walltime::get();
        HOST_PRINT(commShot, "Shot number " << shotNumberEncode << ": Estimate encoded source signals in " << end_t_shot - start_t_shot << " sec.\n");
    }
}

/*! \brief Apply the Wiener filter to a synthetic source
 \param sources Synthetic source
 \param shotInd Shot index of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::applyFilter(KITGPI::Acquisition::Sources<ValueType> &sourcesEncode, IndexType shotNumberEncode, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode) const
{
    IndexType numshotsIncr = sourceSettingsEncode.size();    
    SCAI_ASSERT_ERROR(numshotsIncr != 0, "numshotsIncr == 0 when shotNumberEncode = " + std::to_string(shotNumberEncode));
        
    //get seismogram that corresponds to source type
    lama::DenseMatrix<ValueType> seismo;
    auto sourceType = Acquisition::SeismogramType(sourcesEncode.getSeismogramTypes().getValue(0) - 1);
    seismo = sourcesEncode.getSeismogramHandler().getSeismogram(sourceType).getData();
    lama::DenseMatrix<ComplexValueType> seismoTrans;
    seismoTrans = lama::cast<ComplexValueType>(seismo);
    // scale to power of two
    seismoTrans.resize(seismo.getRowDistributionPtr(), filter.getColDistributionPtr());

    // apply filter in frequency domain
    lama::fft<ComplexValueType>(seismoTrans, 1);

    lama::DenseVector<ComplexValueType> filterTmp;
    lama::DenseVector<ComplexValueType> signalTmp;
    IndexType count = 0;
    for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
        if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumberEncode) {  
            filter.getRow(filterTmp, shotInd);
            SCAI_ASSERT_ERROR(filterTmp.l2Norm() != 0, "filterTmp.l2Norm() == 0 when shotInd = " + std::to_string(shotInd));
            seismoTrans.getRow(signalTmp, count);
            signalTmp *= filterTmp;
            seismoTrans.setRow(signalTmp, count, common::BinaryOp::COPY);
            count++;
        }
    }

    lama::ifft<ComplexValueType>(seismoTrans, 1);
    seismoTrans *= 1.0 / nFFT;

    // return to time domain
    seismoTrans.resize(seismo.getRowDistributionPtr(), seismo.getColDistributionPtr());
    seismo = lama::real(seismoTrans);
    
    sourcesEncode.getSeismogramHandler().getSeismogram(sourceType).getData() = seismo;
}

/*! \brief Correlate the rows of two matrices.
 \param prod Result
 \param A First matrix
 \param B Second matrix
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::matCorr(scai::lama::DenseVector<ComplexValueType> &prod, lama::DenseMatrix<ValueType> const &A, lama::DenseMatrix<ValueType> const &B, scai::IndexType shotInd)
{
    lama::DenseMatrix<ComplexValueType> ATmp;
    lama::DenseMatrix<ComplexValueType> BTmp;

    ATmp = lama::cast<ComplexValueType>(A);
    BTmp = lama::cast<ComplexValueType>(B);

    auto colDist = std::make_shared<dmemo::NoDistribution>(nFFT);
    ATmp.resize(A.getRowDistributionPtr(), colDist);
    BTmp.resize(A.getRowDistributionPtr(), colDist);

    lama::fft<ComplexValueType>(ATmp, 1);
    lama::fft<ComplexValueType>(BTmp, 1);
    ATmp.conj();
    ATmp.binaryOp(ATmp, common::BinaryOp::MULT, BTmp);

    if (useOffsetMutes) {
        ATmp.redistribute(mutes[shotInd].getDistributionPtr(), colDist);
        ATmp.scaleRows(lama::eval<lama::DenseVector<ComplexValueType>>(lama::cast<ComplexValueType>(mutes[shotInd])));
    }

    ATmp.reduce(prod, 1, common::BinaryOp::ADD, common::UnaryOp::COPY);
}

/*! \brief Correlate and sum all components of two receivers.
 \param sum Result
 \param receiversA First receivers
 \param receiversB Second receivers
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::addComponents(lama::DenseVector<ComplexValueType> &sum, KITGPI::Acquisition::Receivers<ValueType> const &receiversA, KITGPI::Acquisition::Receivers<ValueType> const &receiversB, scai::IndexType shotInd, scai::IndexType shotNumber)
{
    lama::DenseVector<ComplexValueType> filterTmp;

    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (receiversA.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            lama::DenseMatrix<ValueType> seismoA;
            lama::DenseMatrix<ValueType> seismoB;
            seismoA = receiversA.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            seismoB = receiversB.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            
            if (readTaper) {
                lama::DenseMatrix<ValueType> seismoA_tmp(seismoA); // tmps needed so real seismograms don't get tapered
                lama::DenseMatrix<ValueType> seismoB_tmp(seismoB);
                
                Taper::Taper2D<ValueType> seismoTaper;
                seismoTaper.init(seismoA.getRowDistributionPtr(), seismoA.getColDistributionPtr(), seismoA.getContextPtr());
                std::string taperNameTmp = taperName + ".shot_" + std::to_string(shotNumber) + ".mtx";
                seismoTaper.read(taperNameTmp);
                seismoTaper.apply(seismoA_tmp);
                seismoTaper.apply(seismoB_tmp);
                
                matCorr(filterTmp, seismoA_tmp, seismoB_tmp, shotInd);
            } else
                matCorr(filterTmp, seismoA, seismoB, shotInd);
            sum += filterTmp;
        }
    }
}

/*! \brief Sets mutes vector needed to exclude offsets greater than a threshold from the source estimation
 \param sources Sources
 \param receivers Receivers
 \param maxOffset Offset threshold
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::calcOffsetMutes(KITGPI::Acquisition::Sources<ValueType> const &sources, KITGPI::Acquisition::Receivers<ValueType> &receivers, ValueType minOffset, ValueType maxOffset, scai::IndexType shotInd, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    lama::DenseVector<ValueType> offsetSign;
    lama::DenseVector<ValueType> offsetTemp;
    lama::DenseVector<ValueType> offset;
    lama::DenseVector<IndexType> sourceIndexVec;
    lama::DenseVector<IndexType> receiverIndexVec;
    IndexType sourceIndex(0);

    // get index of source position
    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (sources.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            sourceIndexVec = sources.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).get1DCoordinates();
            
            SCAI_ASSERT_ERROR(sourceIndexVec.size() == 1, "offset not defined for more than one source position!");

            sourceIndex = sourceIndexVec.getValue(0);
        }
    }

    // get indices of receiver position and calc mutes
    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (receivers.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            receiverIndexVec = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).get1DCoordinates();
            
            Common::calcOffsets(offset, sourceIndex, receiverIndexVec, modelCoordinates);
            
            offsetSign = offset;
            offsetSign.unaryOp(offsetSign, common::UnaryOp::ABS);
            offsetTemp = offsetSign;
            if (maxOffset != 0) {
                offsetSign /= -maxOffset;
                offsetSign.unaryOp(offsetSign, common::UnaryOp::CEIL);
                offsetSign.unaryOp(offsetSign, common::UnaryOp::SIGN);
                offsetSign += 1.0; // offsetSign = 1,1,1,...,1,0,...,0,0,0
            } else {
                offsetSign = 1.0;
            }
            if (minOffset != 0) {
                offsetTemp /= -minOffset;
                offsetTemp.unaryOp(offsetTemp, common::UnaryOp::CEIL);
                offsetTemp.unaryOp(offsetTemp, common::UnaryOp::SIGN);
                offsetTemp *= -1.0; // offsetSign = 0,0,0,...,0,1,...,1,1,1
                offsetSign *= offsetTemp;
            }
            offsets[shotInd] = offset;
            mutes[shotInd] = offsetSign;
            if (useSourceEncode == 0) {
                receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).setOffsets(offsets);
            }
            break;
        }
    }
}

/*! \brief Sets mutes vector needed to exclude offsets greater than a threshold from the source estimation
 \param sources Sources
 \param receivers Receivers
 \param maxOffset Offset threshold
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::calcRefTraces(Configuration::Configuration const &config, scai::IndexType shotInd, KITGPI::Acquisition::Receivers<ValueType> &receivers, Taper::Taper1D<ValueType> const &sourceSignalTaper)
{
    lama::DenseMatrix<ValueType> data;
    ValueType mainVelocity = config.getAndCatch("mainVelocity", 0.0);
    ValueType DT = config.get<ValueType>("DT");
    SCAI_ASSERT_ERROR(mainVelocity != 0.0, "mainVelocity cannot be zero when calculating the referenced trace!");
        
    // get indices of receiver position and calc mutes
    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (receivers.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            data = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            
            IndexType NR = data.getNumRows();
            IndexType NT = data.getNumColumns();
            lama::DenseVector<ValueType> rowTemp(NT, 0.0);
            lama::DenseVector<ValueType> refTrace(NT, 0.0);
            lama::DenseVector<ValueType> refTraceSum(NT, 0.0);
            lama::DenseVector<ValueType> offsetSign;
            lama::DenseVector<ValueType> offsetUse;
            
            IndexType count = 0;
            if (useOffsetMutes) {
                offsetSign = mutes[shotInd] * offsets[shotInd];
                offsetSign.unaryOp(offsetSign, common::UnaryOp::SIGN);
                if (offsetSign.sum() <= 0) {
                    offsetSign *= -1;
                }
                offsetSign += 1;
                offsetSign.unaryOp(offsetSign, common::UnaryOp::SIGN);
                offsetUse = mutes[shotInd] * offsets[shotInd];
                offsetUse = offsetSign * offsetUse; // we use only one side signals
                for (IndexType ir = 0; ir < NR; ir++) {
                    if (offsetUse.getValue(ir) != 0) {
                        count++;
                        data.getRow(rowTemp, ir);
                        IndexType tCut = static_cast<IndexType>(std::abs(offsetUse.getValue(ir)) / mainVelocity / DT + 0.5);
                        for (IndexType it = tCut; it < NT; it++) {
                            refTrace.setValue(it-tCut, rowTemp.getValue(it));
                        }
                        refTraceSum += refTrace;
                    }
                }
            } else {
                for (IndexType ir = 0; ir < NR; ir++) {
                    count++;
                    data.getRow(refTrace, ir);
                    refTraceSum += refTrace;
                }
            }
            SCAI_ASSERT_ERROR(count != 0, "No trace selected for calculating the referenced trace!");
            if (count != 0)
                refTraceSum *= 1.0 / count;
            if (config.get<IndexType>("useSourceSignalTaper") == 2)
                sourceSignalTaper.apply(refTraceSum);
            
            refTraces[iComponent].setRow(refTraceSum, shotInd, common::BinaryOp::COPY);
            receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getRefTraces() = refTraces[iComponent];
            if (useSourceEncode == 0) {
                receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).setOffsets(offsets);
            }
        }
    }
}

/*! \brief Sets mutes vector needed to exclude offsets greater than a threshold from the source estimation
 \param sources Sources
 \param receivers Receivers
 \param maxOffset Offset threshold
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::applyOffsetMute(Configuration::Configuration const &config, scai::IndexType shotInd, KITGPI::Acquisition::Receivers<ValueType> &receivers)
{
    if (useOffsetMutes) {
        for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
            if (receivers.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
                lama::DenseMatrix<ValueType> &data = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
                
                data.scaleRows(mutes[shotInd]);
            }
        }
    }
}

/*! \brief Calculate the Wiener filter
 \param receivers Synthetic receivers
 \param receiversTrue Observed receivers
 \param shotNumber Shot number of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::calcOffsetMutesEncode(scai::dmemo::CommunicatorPtr commShot, scai::IndexType shotNumberEncode, KITGPI::Configuration::Configuration const &config, KITGPI::Workflow::Workflow<ValueType> const &workflow, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, KITGPI::Acquisition::Receivers<ValueType> &receiversEncode)
{
    if (useSourceEncode != 0) {
        double start_t_shot, end_t_shot; /* For timing */
        start_t_shot = common::Walltime::get();
                
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings;  
        Acquisition::Sources<ValueType> sources;
        ValueType shotIncr = 0; 
        sources.getAcquisitionSettings(config, shotIncr); // to get numshots
        bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
        if (useStreamConfig) {
            Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
            Acquisition::Coordinates<ValueType> modelCoordinatesBig(configBig);
            std::vector<Acquisition::coordinate3D> cutCoordinates;
            std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;
            sourceSettingsBig = sources.getSourceSettings(); 
            Acquisition::getCutCoord(config, cutCoordinates, sourceSettingsBig, modelCoordinates, modelCoordinatesBig);
            Acquisition::getSettingsPerShot(sourceSettings, sourceSettingsBig, cutCoordinates, modelCoordinates, config.get<IndexType>("BoundaryWidth"));
        } else {
            sourceSettings = sources.getSourceSettings(); 
        }
        std::vector<scai::IndexType> uniqueShotNos;
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
        IndexType numshotsIncr = sourceSettingsEncode.size();
        std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncodeTemp; // to control the useSourceEncode
        KITGPI::Acquisition::Receivers<ValueType> receivers;
                
        for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
            if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumberEncode) {
                IndexType shotNumber = uniqueShotNos[sourceSettingsEncode[shotInd].row];                 
                std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
                Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
                sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);
                
                if (config.get<IndexType>("useReceiversPerShot") != 0) {
                    receivers.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncodeTemp);
                }
                
                calcOffsetMutes(sources, receivers, workflow.getMinOffset(), workflow.getMaxOffset(), shotInd, modelCoordinates);                                              
            }
        }
        for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
            if (receiversEncode.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
                receiversEncode.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).setOffsets(offsets);
            }
        }
        end_t_shot = common::Walltime::get();
        HOST_PRINT(commShot, "Shot number " << shotNumberEncode << ": Calculate encoded offsets in " << end_t_shot - start_t_shot << " sec.\n");
    }
}

/*! \brief Calculate the Wiener filter
 \param receivers Synthetic receivers
 \param receiversTrue Observed receivers
 \param shotNumber Shot number of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::calcRefTracesEncode(scai::dmemo::CommunicatorPtr commShot, scai::IndexType shotNumberEncode, KITGPI::Configuration::Configuration const &config, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, KITGPI::Acquisition::Receivers<ValueType> &receiversEncode, Taper::Taper1D<ValueType> const &sourceSignalTaper)
{
    if (useSourceEncode != 0) {
        double start_t_shot, end_t_shot; /* For timing */
        start_t_shot = common::Walltime::get();
                
        IndexType numshotsIncr = sourceSettingsEncode.size();
        KITGPI::Acquisition::Receivers<ValueType> receivers;
        IndexType countDecode = 0;
        for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
            if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumberEncode) {
                for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
                    if (receiversEncode.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
                        receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData() = receiversEncode.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getDataDecode(countDecode);
                    }
                }
                
                calcRefTraces(config, shotInd, receivers, sourceSignalTaper);
                countDecode++;
            }
        }
        for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
            if (receiversEncode.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
                receiversEncode.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getRefTraces() = refTraces[iComponent];
            }
        }
        end_t_shot = common::Walltime::get();
        HOST_PRINT(commShot, "Shot number " << shotNumberEncode << ": Calculate encoded refTraces in " << end_t_shot - start_t_shot << " sec.\n");
    }
}

/*! \brief Calculate the Wiener filter
 \param receivers Synthetic receivers
 \param receiversTrue Observed receivers
 \param shotNumber Shot number of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::applyOffsetMuteEncode(scai::dmemo::CommunicatorPtr commShot, scai::IndexType shotNumberEncode, KITGPI::Configuration::Configuration const &config, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, KITGPI::Acquisition::Receivers<ValueType> &receiversEncode)
{
    if (useSourceEncode != 0) {
        double start_t_shot, end_t_shot; /* For timing */
        start_t_shot = common::Walltime::get();
                
        IndexType numshotsIncr = sourceSettingsEncode.size();
        KITGPI::Acquisition::Receivers<ValueType> receivers;
        IndexType countDecode = 0;
        for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
            if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumberEncode) {
                for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
                    if (receiversEncode.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
                        receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData() = receiversEncode.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getDataDecode(countDecode);
                        
                        applyOffsetMute(config, shotInd, receivers);
                        
                        receiversEncode.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getDataDecode(countDecode) = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
                    }
                }                
                countDecode++;
            }
        }
        end_t_shot = common::Walltime::get();
        HOST_PRINT(commShot, "Shot number " << shotNumberEncode << ": Apply encoded offsetMute in " << end_t_shot - start_t_shot << " sec.\n");
    }
}

/*! \brief Set the refTrace to a synthetic source
 \param sources Synthetic source
 \param shotInd Shot index of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::setRefTracesToSource(KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> const &receivers, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> sourceSettingsEncode, scai::IndexType shotIndTrue, scai::IndexType shotNumberEncode)
{
    //get seismogram that corresponds to source type
    auto sourceType = Acquisition::SeismogramType(sources.getSeismogramTypes().getValue(0) - 1);
    
    lama::DenseMatrix<ValueType> data;
    lama::DenseMatrix<ValueType> refTracesTemp;
    lama::DenseVector<ValueType> refTrace;
    refTracesTemp = receivers.getSeismogramHandler().getSeismogram(sourceType).getRefTraces();
    data = sources.getSeismogramHandler().getSeismogram(sourceType).getData();
    if (useSourceEncode == 0) {
        refTracesTemp.getRow(refTrace, shotIndTrue);
        data.setRow(refTrace, 0, common::BinaryOp::COPY);
    } else {
        IndexType numshotsIncr = sourceSettingsEncode.size();
        IndexType count = 0;
        for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
            if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumberEncode) {
                refTracesTemp.getRow(refTrace, shotInd);
                data.setRow(refTrace, count, common::BinaryOp::COPY);
                count++;
            }
        }
    }
    sources.getSeismogramHandler().getSeismogram(sourceType).getData() = data;
}

template class KITGPI::SourceEstimation<double>;
template class KITGPI::SourceEstimation<float>;
