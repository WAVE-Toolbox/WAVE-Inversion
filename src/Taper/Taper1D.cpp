#include "Taper1D.hpp"
#include <IO/IO.hpp>

using namespace scai;

/*! \brief Initialize taper
 \param dist Distribution
 \param ctx Context
 \param dir Direction (0=taper columns, 1=taper rows)
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, bool dir)
{
    direction = dir;
    data.allocate(dist);
    data = 1.0; // in this state the taper does nothing when applied
    data.setContextPtr(ctx);
}

/*! \brief Get direction of taper
 */
template <typename ValueType>
bool KITGPI::Taper::Taper1D<ValueType>::getDirection() const
{
    return (direction);
}

/*! \brief to calculate a time damping taper for seismogram
 \param timeDampingFactor timeDampingFactor
 \param DT DT
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::calcTimeDampingTaper(ValueType timeDampingFactor, ValueType DT)
{
    data = lama::linearDenseVector<ValueType>(data.size(), 0.0, -DT);
    data *= timeDampingFactor;
    data.unaryOp(data, common::UnaryOp::EXP);
}

/*! \brief Wrapper to calculate a cosine taper with one transition zone
 \param iStart Start index of transition zone
 \param iEnd End index of transition zone
 \param reverse 0 = Taper starts at 0 and ends with 1 (slope >= 0), 1 = Taper starts at 1 and ends with 0 (slope <= 0)
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::calcCosineTaper(IndexType iStart, IndexType iEnd, bool reverse)
{
    SCAI_ASSERT_ERROR(iStart >= 0 && iEnd < data.size() && iStart < iEnd, "invalid taper edges");

    if (reverse)
        calcCosineTaperDown(data, iStart, iEnd);
    else
        calcCosineTaperUp(data, iStart, iEnd);
}

/*! \brief Wrapper to calculate a cosine taper with two transition zones
 \param iStart1 Start index of first transition zone
 \param iEnd1 End index of first transition zone
 \param iStart2 Start index of second transition zone
 \param iEnd2 End index of second transition zone
 \param reverse 0 = Taper starts at 0 and ends with 0 , 1 = Taper starts at 1 and ends with 1
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::calcCosineTaper(IndexType iStart1, IndexType iEnd1, IndexType iStart2, IndexType iEnd2, bool reverse)
{
    SCAI_ASSERT_ERROR(iStart1 >= 0 && iEnd2 < data.size() && iStart1 < iEnd1 && iStart2 < iEnd2 && iEnd1 <= iStart2, "invalid taper edges");

    lama::DenseVector<ValueType> helpTaper;
    calcCosineTaperUp(helpTaper, iStart1, iEnd1);
    calcCosineTaperDown(data, iStart2, iEnd2);
    data *= helpTaper;

    if (reverse)
        data = -data + 1;
}

/*! \brief Wrapper to calculate a cosine taper using the inverted source
 \param seismograms Seismogram handler of the inverted source
 \param lowerCornerFreq Lower corner frequency
 \param upperCornerFreq Upper corner frequency
 \param DT DT
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::calcCosineTaper(KITGPI::Acquisition::SeismogramHandler<ValueType> const &seismograms, ValueType lowerCornerFreq, ValueType upperCornerFreq, ValueType DT, scai::hmemo::ContextPtr ctx)
{
    bool isSeismic = seismograms.getIsSeismic();
    Acquisition::Seismogram<ValueType> thisSeismogram;
    for (scai::IndexType iComponent = 0; iComponent < KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (isSeismic && seismograms.getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            thisSeismogram = seismograms.getSeismogram(Acquisition::SeismogramType(iComponent));
        } else if (!isSeismic && seismograms.getNumTracesGlobal(Acquisition::SeismogramTypeEM(iComponent)) != 0) {
            thisSeismogram = seismograms.getSeismogram(Acquisition::SeismogramTypeEM(iComponent));
        }
        if (thisSeismogram.getData().getNumRows()!=0)
            break;
    }
    scai::lama::DenseMatrix<ValueType> tempData = thisSeismogram.getData();
    scai::lama::DenseVector<ValueType> tempRow;
    Common::calcEnvelope(tempData);
    tempData.getRow(tempRow, 0);
    ValueType maxValue = tempRow.max();
    IndexType maxIndex = 0;
    IndexType tStepEnd = tempRow.size();
    for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
        if (tempRow[tStep] == maxValue) {
            maxIndex = tStep;
            break;
        }
    }
    ValueType medFreq = (lowerCornerFreq + upperCornerFreq) / 2;
    IndexType period = floor(1.0 / medFreq / DT);
    IndexType iStart1 = maxIndex - period;
    IndexType iEnd1 = maxIndex - period / 2;
    IndexType iStart2 = maxIndex + period / 2;
    IndexType iEnd2 = maxIndex + period;
    if (iStart1 < 0)
        iStart1 = 0;
    if (iEnd1 < 10)
        iEnd1 = 10;
    SCAI_ASSERT_ERROR(iEnd2 < tStepEnd, "iEnd2 >= tStepEnd");
    this->init(std::make_shared<dmemo::NoDistribution>(tStepEnd), ctx, 1);
    this->calcCosineTaper(iStart1, iEnd1, iStart2, iEnd2, 0);
}

/*! \brief Calculate cosine taper which starts with 0 and ends with 1 (slope >= 0)
 * \param result Result vector
 \param iStart Start index of transition zone
 \param iEnd End index of transition zone
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::calcCosineTaperUp(lama::DenseVector<ValueType> &result, IndexType iStart, IndexType iEnd)
{
    // first part of taper
    lama::DenseVector<ValueType> firstPart(iStart, 0.0);
    lama::DenseVector<ValueType> tmpResult;

    // second part of taper
    lama::DenseVector<ValueType> secondPart = lama::linearDenseVector<ValueType>(iEnd - iStart, 1.0, -1.0 / (iEnd - iStart));
    secondPart *= M_PI / 2.0;
    secondPart.unaryOp(secondPart, common::UnaryOp::COS);
    secondPart.binaryOpScalar(secondPart, 2.0, common::BinaryOp::POW, false);
    tmpResult.cat(firstPart, secondPart);

    // third part of taper
    firstPart = tmpResult;
    secondPart.allocate(data.size() - iEnd);
    secondPart = 1.0;
    tmpResult.cat(firstPart, secondPart);

    result.assignDistribute(tmpResult, data.getDistributionPtr());
}

/*! \brief Calculate cosine taper which starts with 1 and ends with 0 (slope <= 0)
 \param result Result vector
 \param iStart Start index of transition zone
 \param iEnd End index of transition zone
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::calcCosineTaperDown(lama::DenseVector<ValueType> &result, IndexType iStart, IndexType iEnd)
{
    // first part of taper
    lama::DenseVector<ValueType> firstPart(iStart, 1.0);
    lama::DenseVector<ValueType> tmpResult;

    // second part of taper
    lama::DenseVector<ValueType> secondPart = lama::linearDenseVector<ValueType>(iEnd - iStart, 0.0, 1.0 / (iEnd - iStart));
    secondPart *= M_PI / 2.0;
    secondPart.unaryOp(secondPart, common::UnaryOp::COS);
    secondPart.binaryOpScalar(secondPart, 2.0, common::BinaryOp::POW, false);
    tmpResult.cat(firstPart, secondPart);

    // third part of taper
    firstPart = tmpResult;
    secondPart.allocate(data.size() - iEnd);
    secondPart = 0.0;
    tmpResult.cat(firstPart, secondPart);

    result.assignDistribute(tmpResult, data.getDistributionPtr());
}

/*! \brief Wrapper to support SeismogramHandler
 \param seismograms SeismogramHandler object
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::apply(KITGPI::Acquisition::SeismogramHandler<ValueType> &seismograms) const
{
    bool isSeismic = seismograms.getIsSeismic();
    for (scai::IndexType iComponent = 0; iComponent < KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        if (isSeismic && seismograms.getNumTracesGlobal(Acquisition::SeismogramType(iComponent)) != 0) {
            Acquisition::Seismogram<ValueType> &thisSeismogram = seismograms.getSeismogram(Acquisition::SeismogramType(iComponent));
            apply(thisSeismogram);
        } else if (!isSeismic && seismograms.getNumTracesGlobal(Acquisition::SeismogramTypeEM(iComponent)) != 0) {
            Acquisition::Seismogram<ValueType> &thisSeismogram = seismograms.getSeismogram(Acquisition::SeismogramTypeEM(iComponent));
            apply(thisSeismogram);
        }
    }
}

/*! \brief Apply taper to a single seismogram
 \param seismogram Seismogram
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::apply(KITGPI::Acquisition::Seismogram<ValueType> &seismogram) const
{
    lama::DenseMatrix<ValueType> &seismogramData = seismogram.getData();
    apply(seismogramData);
}

/*! \brief Apply taper to a DenseMatrix
 \param mat Seismogram
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::apply(lama::DenseMatrix<ValueType> &mat) const
{
    if (direction == 0)
        mat.scaleRows(data); // scaleRows means, that every row is scaled with one entry in data
    else
        mat.scaleColumns(data);
}

/*! \brief Apply taper to a Gradient
 \param grad Gradient
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::apply(scai::lama::DenseVector<ValueType> &trace) const
{
    trace *= data;
}

/*! \brief Apply taper to a Gradient
 \param grad Gradient
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::apply(KITGPI::Gradient::Gradient<ValueType> &grad) const
{
    grad *= data;
}

/*! \brief Read a taper from file
 */
template <typename ValueType>
void KITGPI::Taper::Taper1D<ValueType>::read(std::string filename, IndexType fileFormat)
{
    IO::readVector(data, filename, fileFormat);
}

template class KITGPI::Taper::Taper1D<double>;
template class KITGPI::Taper::Taper1D<float>;
