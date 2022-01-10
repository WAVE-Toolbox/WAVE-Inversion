#include "SourceEstimation.hpp"

using namespace scai;

/*! \brief Initialize based on config which also initializes sourceSignalTaper
 \param config Configuration
 \param ctx Context Ptr
 \param sourceDistribution Source distribution pointer
 \param sourceSignalTaper 1D TaperPtr
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr sourceDistribution, Taper::Taper1D<ValueType> &sourceSignalTaper)
{
    std::string equationType = config.get<std::string>("equationType"); 
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);  
    isSeismic = Common::checkEquationType<ValueType>(equationType);
    
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    if (config.get<IndexType>("useSeismogramTaper") != 0 && config.get<IndexType>("useSeismogramTaper") != 2) {
        if (config.get<IndexType>("useSeismogramTaper") == 4) {        
            init(tStepEnd, sourceDistribution, config.get<ValueType>("waterLevel"), config.get<std::string>("seismogramTaperName") + ".SrcEst");
        } else {     
            init(tStepEnd, sourceDistribution, config.get<ValueType>("waterLevel"), config.get<std::string>("seismogramTaperName"));
        }
    } else {
        init(tStepEnd, sourceDistribution, config.get<ValueType>("waterLevel"));
    }
    
    if (config.getAndCatch("minOffsetSrcEst", 0.0) || config.get<bool>("maxOffsetSrcEst"))
        useOffsetMutes = true;

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
//     waterLevel = common::Math::pow<ValueType>(waterLvl, 2.0) * nFFT;
    waterLevel = waterLvl;
    filter.allocate(sourceDistribution, std::make_shared<dmemo::NoDistribution>(nFFT));
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
void KITGPI::SourceEstimation<ValueType>::estimateSourceSignal(KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Receivers<ValueType> &receiversTrue, IndexType shotInd, IndexType shotNr)
{
    lama::DenseVector<ComplexValueType> filterTmp1(filter.getColDistributionPtr(), 0.0);
    lama::DenseVector<ComplexValueType> filterTmp2(filter.getColDistributionPtr(), 0.0);

    addComponents(filterTmp1, receivers, receivers, shotNr);
//     filterTmp1 += waterLevel * receivers.getSeismogramHandler().getNumTracesTotal();
    filterTmp1 += waterLevel * filterTmp1.maxNorm();
    filterTmp1.unaryOp(filterTmp1, common::UnaryOp::RECIPROCAL);
    addComponents(filterTmp2, receivers, receiversTrue, shotNr);
    filterTmp1 *= filterTmp2;
    
    filter.setRow(filterTmp1, shotInd, common::BinaryOp::COPY);
}

/*! \brief Apply the Wiener filter to a synthetic source
 \param sources Synthetic source
 \param shotInd Shot index of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::applyFilter(KITGPI::Acquisition::Sources<ValueType> &sources, IndexType shotInd) const
{
    //get seismogram that corresponds to source type
    lama::DenseMatrix<ValueType> seismo;
    if (isSeismic) {
        auto sourceType = Acquisition::SeismogramType(sources.getSeismogramTypes().getValue(0) - 1);
        seismo = sources.getSeismogramHandler().getSeismogram(sourceType).getData();
    } else {
        auto sourceType = Acquisition::SeismogramTypeEM(sources.getSeismogramTypes().getValue(0) - 1);
        seismo = sources.getSeismogramHandler().getSeismogram(sourceType).getData();
    }
    lama::DenseMatrix<ComplexValueType> seismoTrans;
    seismoTrans = lama::cast<ComplexValueType>(seismo);

    // scale to power of two
    seismoTrans.resize(seismo.getRowDistributionPtr(), filter.getColDistributionPtr());

    // apply filter in frequency domain
    lama::fft<ComplexValueType>(seismoTrans, 1);

    lama::DenseVector<ComplexValueType> filterTmp;
    filter.getRow(filterTmp, shotInd);
    seismoTrans.scaleColumns(filterTmp);

    lama::ifft<ComplexValueType>(seismoTrans, 1);
    seismoTrans *= 1.0 / nFFT;

    // return to time domain
    seismoTrans.resize(seismo.getRowDistributionPtr(), seismo.getColDistributionPtr());
    seismo = lama::real(seismoTrans);
    
    if (isSeismic) {
        auto sourceType = Acquisition::SeismogramType(sources.getSeismogramTypes().getValue(0) - 1);
        sources.getSeismogramHandler().getSeismogram(sourceType).getData() = seismo;
    } else {
        auto sourceType = Acquisition::SeismogramTypeEM(sources.getSeismogramTypes().getValue(0) - 1);
        sources.getSeismogramHandler().getSeismogram(sourceType).getData() = seismo;
    }
}

/*! \brief Correlate the rows of two matrices.
 \param prod Result
 \param A First matrix
 \param B Second matrix
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::matCorr(lama::DenseVector<ComplexValueType> &prod, lama::DenseMatrix<ValueType> const &A, lama::DenseMatrix<ValueType> const &B, IndexType iComponent)
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
        ATmp.scaleRows(lama::eval<lama::DenseVector<ComplexValueType>>(lama::cast<ComplexValueType>(mutes[iComponent])));
    }

    ATmp.reduce(prod, 1, common::BinaryOp::ADD, common::UnaryOp::COPY);
}

/*! \brief Correlate and sum all components of two receivers.
 \param sum Result
 \param receiversA First receivers
 \param receiversB Second receivers
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::addComponents(lama::DenseVector<ComplexValueType> &sum, KITGPI::Acquisition::Receivers<ValueType> const &receiversA, KITGPI::Acquisition::Receivers<ValueType> const &receiversB, IndexType shotNumber)
{
    lama::DenseVector<ComplexValueType> filterTmp;

    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (receiversA.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            lama::DenseMatrix<ValueType> seismoA;
            lama::DenseMatrix<ValueType> seismoB;
            if (isSeismic) {
                seismoA = receiversA.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
                seismoB = receiversB.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            } else {
                seismoA = receiversA.getSeismogramHandler().getSeismogram(Acquisition::SeismogramTypeEM(iComponent)).getData();
                seismoB = receiversB.getSeismogramHandler().getSeismogram(Acquisition::SeismogramTypeEM(iComponent)).getData();
            }
            if (readTaper) {
                lama::DenseMatrix<ValueType> seismoA_tmp(seismoA); // tmps needed so real seismograms don't get tapered
                lama::DenseMatrix<ValueType> seismoB_tmp(seismoB);
                
                Taper::Taper2D<ValueType> seismoTaper;
                seismoTaper.init(seismoA.getRowDistributionPtr(), seismoA.getColDistributionPtr(), seismoA.getContextPtr());
                std::string taperNameTmp = taperName + ".shot_" + std::to_string(shotNumber) + ".mtx";
                seismoTaper.read(taperNameTmp);
                seismoTaper.apply(seismoA_tmp);
                seismoTaper.apply(seismoB_tmp);
                
                matCorr(filterTmp, seismoA_tmp, seismoB_tmp, iComponent);
            } else
                matCorr(filterTmp, seismoA, seismoB, iComponent);
            sum += filterTmp;
        }
    }
}

/*! \brief Sets mutes vector needed to exclude offsets greater than a threshold from the source estimation
 \param sources Sources
 \param receivers Receivers
 \param maxOffset Offset threshold
 \param NX Number of grid points in x-direction
 \param NY Number of grid points in y-direction
 \param NZ Number of grid points in z-direction
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::calcOffsetMutes(KITGPI::Acquisition::Sources<ValueType> const &sources, KITGPI::Acquisition::Receivers<ValueType> &receivers, ValueType minOffset, ValueType maxOffset, KITGPI::Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    lama::DenseVector<ValueType> offsetsSign;
    lama::DenseVector<ValueType> offsetsTemp;
    lama::DenseVector<IndexType> sourceIndexVec;
    lama::DenseVector<IndexType> receiverIndexVec;
    IndexType sourceIndex(0);

    bool sourceIndexSet(false);
    bool isSeismic = sources.getSeismogramHandler().getIsSeismic();

    // get index of source position
    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (sources.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramTypeEM(iComponent)) != 0) {
            if (isSeismic)
                sourceIndexVec = sources.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).get1DCoordinates();
            else
                sourceIndexVec = sources.getSeismogramHandler().getSeismogram(Acquisition::SeismogramTypeEM(iComponent)).get1DCoordinates();

            if (sourceIndexVec.size() > 1 || (sourceIndexSet && sourceIndexVec.getValue(0) != sourceIndex))
                COMMON_THROWEXCEPTION("offset not defined for more than one source position");

            sourceIndex = sourceIndexVec.getValue(0);
        }
    }

    // get indices of receiver position and calc mutes
    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (receivers.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramTypeEM(iComponent)) != 0) {
            if (isSeismic)
                receiverIndexVec = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).get1DCoordinates();
            else
                receiverIndexVec = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramTypeEM(iComponent)).get1DCoordinates();
            Common::calcOffsets(offsets, sourceIndex, receiverIndexVec, modelCoordinates);
            
            offsetsSign = offsets;
            offsetsSign.unaryOp(offsetsSign, common::UnaryOp::ABS);
            offsetsTemp = offsetsSign;
            if (maxOffset != 0) {
                offsetsSign /= -maxOffset;
                offsetsSign.unaryOp(offsetsSign, common::UnaryOp::CEIL);
                offsetsSign.unaryOp(offsetsSign, common::UnaryOp::SIGN);
                offsetsSign += 1.0; // offsetsSign = 1,1,1,...,1,0,...,0,0,0
            } else {
                offsetsSign = 1.0;
            }
            if (minOffset != 0) {
                offsetsTemp /= -minOffset;
                offsetsTemp.unaryOp(offsetsTemp, common::UnaryOp::CEIL);
                offsetsTemp.unaryOp(offsetsTemp, common::UnaryOp::SIGN);
                offsetsTemp *= -1.0; // offsetsSign = 0,0,0,...,0,1,...,1,1,1
                offsetsSign *= offsetsTemp;
            }

            mutes[iComponent] = offsetsSign;
            if (isSeismic)
                receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).setOffset(offsets);
            else
                receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramTypeEM(iComponent)).setOffset(offsets);
        }
    }
}

/*! \brief Sets mutes vector needed to exclude offsets greater than a threshold from the source estimation
 \param sources Sources
 \param receivers Receivers
 \param maxOffset Offset threshold
 \param NX Number of grid points in x-direction
 \param NY Number of grid points in y-direction
 \param NZ Number of grid points in z-direction
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::calcRefTrace(Configuration::Configuration const &config, KITGPI::Acquisition::Receivers<ValueType> &receivers, Taper::Taper1D<ValueType> const &sourceSignalTaper)
{
    lama::DenseMatrix<ValueType> data;
    ValueType mainVelocity = config.getAndCatch("mainVelocity", 0.0);
    ValueType DT = config.get<ValueType>("DT");
    SCAI_ASSERT_ERROR(mainVelocity != 0.0, "mainVelocity cannot be zero when calculating the referenced trace!");
    
    bool isSeismic = receivers.getSeismogramHandler().getIsSeismic();

    // get indices of receiver position and calc mutes
    for (IndexType iComponent = 0; iComponent < Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (receivers.getSeismogramHandler().getNumTracesGlobal(Acquisition::SeismogramTypeEM(iComponent)) != 0) {
            if (isSeismic)
                data = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).getData();
            else
                data = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramTypeEM(iComponent)).getData();
            
            IndexType NX = data.getNumRows();
            IndexType NT = data.getNumColumns();
            lama::DenseVector<ValueType> rowTemp(NT, 0.0);
            lama::DenseVector<ValueType> refTrace(NT, 0.0);
            lama::DenseVector<ValueType> refTraceSum(NT, 0.0);
            lama::DenseVector<ValueType> refTrace1;
            lama::DenseVector<ValueType> offsetSign;
            lama::DenseVector<ValueType> offsetUse;
            IndexType icount = 0;
            offsetSign = mutes[iComponent] * offsets;
            offsetSign.unaryOp(offsetSign, common::UnaryOp::SIGN);
            if (offsetSign.sum() <= 0) {
                offsetSign *= -1;
            }
            offsetSign += 1;
            offsetSign.unaryOp(offsetSign, common::UnaryOp::SIGN);
            offsetUse = mutes[iComponent] * offsets;
            offsetUse = offsetSign * offsetUse; // we use only one side signals
            
            for (IndexType ix = 0; ix < NX; ix++) {
                if (offsetUse.getValue(ix) != 0) {
                    icount++;
                    data.getRow(rowTemp, ix);
                    IndexType tCut = static_cast<IndexType>(std::abs(offsetUse.getValue(ix)) / mainVelocity / DT + 0.5);
                    for (IndexType it = tCut; it < NT; it++) {
                        refTrace.setValue(it-tCut, rowTemp.getValue(it));
                    }
                    if (icount == 1) {
                        refTrace1 = refTrace;
                    } else {
                        refTrace *= refTrace1.maxNorm() / refTrace.maxNorm();
                    }
                    refTraceSum += refTrace;
                }
            }
            SCAI_ASSERT_ERROR(icount != 0, "No trace selected for calculating the referenced trace!");
            if (icount != 0)
                refTraceSum *= 1.0 / icount;
            if (config.get<IndexType>("useSourceSignalTaper") == 2)
                sourceSignalTaper.apply(refTraceSum);
            
            if (isSeismic) {
                receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).setRefTrace(refTraceSum);
                receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType(iComponent)).setOffset(offsets);
            } else {
                receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramTypeEM(iComponent)).setRefTrace(refTraceSum);
                receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramTypeEM(iComponent)).setOffset(offsets);
            }
        }
    }
}

/*! \brief Set the refTrace to a synthetic source
 \param sources Synthetic source
 \param shotInd Shot index of source
 */
template <typename ValueType>
void KITGPI::SourceEstimation<ValueType>::setRefTraceToSource(KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Acquisition::Receivers<ValueType> const &receivers)
{
    //get seismogram that corresponds to source type
    lama::DenseMatrix<ValueType> seismo;
    lama::DenseVector<ValueType> refTrace;
    if (isSeismic) {
        auto sourceType = Acquisition::SeismogramType(sources.getSeismogramTypes().getValue(0) - 1);
        seismo = sources.getSeismogramHandler().getSeismogram(sourceType).getData();
        refTrace = receivers.getSeismogramHandler().getSeismogram(sourceType).getRefTrace();
    } else {
        auto sourceType = Acquisition::SeismogramTypeEM(sources.getSeismogramTypes().getValue(0) - 1);
        seismo = sources.getSeismogramHandler().getSeismogram(sourceType).getData();
        refTrace = receivers.getSeismogramHandler().getSeismogram(sourceType).getRefTrace();
    }
    for (IndexType i = 0; i < seismo.getNumRows(); i++) {
        seismo.setRow(refTrace, i, common::BinaryOp::COPY);
    }
    if (isSeismic) {
        auto sourceType = Acquisition::SeismogramType(sources.getSeismogramTypes().getValue(0) - 1);
        sources.getSeismogramHandler().getSeismogram(sourceType).getData() = seismo;
    } else {
        auto sourceType = Acquisition::SeismogramTypeEM(sources.getSeismogramTypes().getValue(0) - 1);
        sources.getSeismogramHandler().getSeismogram(sourceType).getData() = seismo;
    }
}

template class KITGPI::SourceEstimation<double>;
template class KITGPI::SourceEstimation<float>;
