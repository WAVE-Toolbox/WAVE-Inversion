#include "FK.hpp"
using namespace scai;

/*! \brief Initialize the filter transfer function. In this state nothing is filtered when applying it.
 \param dt Sampling period of the signal the filter should be applied on
 \param nt Length of the signal the filter should be applied on
 */
template <typename ValueType>
void KITGPI::FK<ValueType>::init(ValueType dt, scai::IndexType nt, ValueType fc, ValueType vmin)
{
    SCAI_ASSERT_ERROR(dt != 0.0, "Can't initialize fk with dt = 0.0")
    SCAI_ASSERT_ERROR(nt != 0, "Can't initialize fk with nt = 0")
    SCAI_ASSERT_ERROR(fc != 0.0, "Can't initialize fk with fc = 0.0")
    SCAI_ASSERT_ERROR(vmin != 0.0, "Can't initialize fk with vmin = 0.0")
    scai::IndexType len = Common::calcNextPowTwo<ValueType>(nt - 1);
    ValueType df = 1.0 / (len * dt);
    ValueType fNyquist = 1.0 / (2.0 * dt);
    long nFreq = fNyquist / df;
    scai::lama::DenseVector<ValueType> fPos = scai::lama::linearDenseVector<ValueType>(nFreq + 1, 0.0, df);
    scai::lama::DenseVector<ValueType> fNeg = scai::lama::linearDenseVector<ValueType>(nFreq - 1, -(nFreq - 1) * df, df);
    freqVec.cat(fPos, fNeg);
    fc2 = fc * 4.0;
    ValueType kmax = fc2 / vmin;
    ValueType dk = 2.0 * kmax / ValueType(NK-1);
    kVec = scai::lama::linearDenseVector<ValueType>(NK, -kmax, dk);
}

/*! \brief Calculate the transfer function of a FK transform.
 \param filterType Type of filter: "lp" = low pass, "hp" = high pass
 \param fc1 Lower corner frequency in Hz
 \param fc2 Upper corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::FK<ValueType>::calcFKOperatorL(scai::lama::DenseVector<ValueType> offset, scai::lama::DenseMatrix<ComplexValueType> &L) const
{ 
    IndexType NX = offset.size();
    IndexType NK = kVec.size();
    auto distNK = std::make_shared<scai::dmemo::NoDistribution>(NK);
    auto distNX = std::make_shared<scai::dmemo::NoDistribution>(NX);
    L.clear();
    L.allocate(distNK, distNX);
    ComplexValueType j(0.0, 1.0); 
    scai::lama::DenseVector<ComplexValueType> temp;
    for (int ix = 0; ix < NX; ix++) {
        temp = scai::lama::cast<ComplexValueType>(kVec);
        temp *= -j * 2.0 * M_PI * offset.getValue(ix);
        temp.unaryOp(temp, common::UnaryOp::EXP);
        L.setColumn(temp, ix, common::BinaryOp::COPY);
    }
}

/*! \brief Calculate the transfer function of a FK transform.
 \param filterType Type of filter: "lp" = low pass, "hp" = high pass
 \param fc1 Lower corner frequency in Hz
 \param fc2 Upper corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::FK<ValueType>::calcFKOperatorLinv(scai::lama::DenseVector<ValueType> offset, scai::lama::DenseMatrix<ComplexValueType> &Linv) const
{ 
    IndexType NX = offset.size();
    IndexType NK = kVec.size();
    auto distNK = std::make_shared<scai::dmemo::NoDistribution>(NK);
    auto distNX = std::make_shared<scai::dmemo::NoDistribution>(NX);
    Linv.clear();
    Linv.allocate(distNX, distNK);
    ComplexValueType j(0.0, 1.0); 
    scai::lama::DenseVector<ComplexValueType> temp;
    for (int ix = 0; ix < NX; ix++) {
        temp = scai::lama::cast<ComplexValueType>(kVec);
        temp *= j * 2.0 * M_PI * offset.getValue(ix);
        temp.unaryOp(temp, common::UnaryOp::EXP);
        Linv.setRow(temp, ix, common::BinaryOp::COPY);
    }
}

/*! \brief Calculate the transfer function of a FK transform.
 \param filterType Type of filter: "lp" = low pass, "hp" = high pass
 \param fc1 Lower corner frequency in Hz
 \param fc2 Upper corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::FK<ValueType>::FKTransform(scai::lama::DenseMatrix<ValueType> const signal, scai::lama::DenseMatrix<ComplexValueType> &fk, scai::lama::DenseVector<ValueType> const offset) const
{
    scai::IndexType len = freqVec.size();

    scai::lama::DenseMatrix<ComplexValueType> fSignal;
    fSignal = scai::lama::cast<ComplexValueType>(signal);
    fSignal.resize(signal.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(len));
    scai::lama::fft<ComplexValueType>(fSignal, 1);
    
    scai::lama::DenseMatrix<ComplexValueType> L;
    calcFKOperatorL(offset, L);
    IndexType NX = offset.size();

    IndexType indexFc1 = 0;
    IndexType indexFc2 = 0;
    for (int jf = 0; jf < len-1; jf++) {
        if (freqVec.getValue(jf)<=fc1 && freqVec.getValue(jf+1)>=fc1) {
            indexFc1 = jf;
        } else if (freqVec.getValue(jf)<=fc2 && freqVec.getValue(jf+1)>=fc2) {
            indexFc2 = jf;
            break;
        }
    }
    
    IndexType NF = indexFc2 - indexFc1 + 1;
    scai::lama::DenseVector<ComplexValueType> k;
    scai::lama::DenseVector<ComplexValueType> x;
    auto distNK = std::make_shared<scai::dmemo::NoDistribution>(NK);
    auto distNX = std::make_shared<scai::dmemo::NoDistribution>(NX);
    fk.clear();
    fk.allocate(distNK, std::make_shared<scai::dmemo::NoDistribution>(NF));
    for (int jf = 0; jf < NF; jf++) {
        fSignal.getColumn(x, indexFc1+jf);
        k = L * x;
        fk.setColumn(k, jf, common::BinaryOp::COPY);
    }
}

/*! \brief Calculate the transfer function of a FK transform.
 \param filterType Type of filter: "lp" = low pass, "hp" = high pass
 \param fc1 Lower corner frequency in Hz
 \param fc2 Upper corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::FK<ValueType>::inverseFKTransform(scai::lama::DenseMatrix<ValueType> &signal, scai::lama::DenseMatrix<ComplexValueType> const fk, scai::lama::DenseVector<ValueType> const offset) const
{
    scai::IndexType len = freqVec.size();

    scai::lama::DenseMatrix<ComplexValueType> fSignal;
    fSignal.allocate(signal.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(len));
    
    scai::lama::DenseMatrix<ComplexValueType> Linv;
    calcFKOperatorLinv(offset, Linv);
    IndexType NX = offset.size();

    IndexType indexFc1 = 0;
    IndexType indexFc2 = 0;
    for (int jf = 0; jf < len-1; jf++) {
        if (freqVec.getValue(jf)<=fc1 && freqVec.getValue(jf+1)>=fc1) {
            indexFc1 = jf;
        } else if (freqVec.getValue(jf)<=fc2 && freqVec.getValue(jf+1)>=fc2) {
            indexFc2 = jf;
            break;
        }
    }
    
    IndexType NF = indexFc2 - indexFc1 + 1;
    scai::lama::DenseVector<ComplexValueType> x;
    scai::lama::DenseVector<ComplexValueType> k;
    auto distNK = std::make_shared<scai::dmemo::NoDistribution>(NK);
    auto distNX = std::make_shared<scai::dmemo::NoDistribution>(NX);
    if (indexFc1 == 0) { // fc1 = 0
        fk.getColumn(k, 0);
        x = Linv * k;
        fSignal.setColumn(x, 0, common::BinaryOp::COPY);
        for (int jf = 1; jf < NF; jf++) {
            fk.getColumn(k, jf);
            x = Linv * k;
            fSignal.setColumn(x, indexFc1+jf, common::BinaryOp::COPY);
            x = scai::lama::conj(x);
            fSignal.setColumn(x, len-indexFc1-jf, common::BinaryOp::COPY);
        }
    } else {
        for (int jf = 0; jf < NF; jf++) {
            fk.getColumn(k, jf);
            x = Linv * k;
            fSignal.setColumn(x, indexFc1+jf, common::BinaryOp::COPY);
            x = scai::lama::conj(x);
            fSignal.setColumn(x, len-indexFc1-jf, common::BinaryOp::COPY);
        }
    }
    
    fSignal *= (1.0 / ValueType(NK)); // proper fft normalization
    fSignal *= (1.0 / ValueType(len)); // proper fft normalization

    scai::lama::ifft<ComplexValueType>(fSignal, 1);
    fSignal.resize(signal.getRowDistributionPtr(), signal.getColDistributionPtr());
    signal = scai::lama::real(fSignal);
}

template class KITGPI::FK<double>;
template class KITGPI::FK<float>;
