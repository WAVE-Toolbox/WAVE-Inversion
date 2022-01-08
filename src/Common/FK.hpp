#pragma once
#include <scai/common/Complex.hpp>
#include <scai/common/UnaryOp.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/fft.hpp>

#include <algorithm>
#include <string>

#include <Common/Common.hpp>
#include <complex>

namespace KITGPI
{
    //! \brief Class to handle frequency filtering.
    template <typename ValueType>
    class FK
    {
        public:
            //! Default constructor
            FK(){};
            //! Default destructor
            ~FK(){};

            typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;
            
            void init(ValueType dt, scai::IndexType nt, ValueType fmax, ValueType vmin);

            void FKTransform(scai::lama::DenseMatrix<ValueType> const signal, scai::lama::DenseMatrix<ComplexValueType> &fk, scai::lama::DenseVector<ValueType> const offset) const;
            void inverseFKTransform(scai::lama::DenseMatrix<ValueType> &signal, scai::lama::DenseMatrix<ComplexValueType> const fk, scai::lama::DenseVector<ValueType> const offset) const;
        
        private:
            void calcFKOperatorL(scai::lama::DenseVector<ValueType> offset, scai::lama::DenseMatrix<ComplexValueType> &L) const;
            void calcFKOperatorLinv(scai::lama::DenseVector<ValueType> offset, scai::lama::DenseMatrix<ComplexValueType> &Linv) const;
            
            ValueType fc1 = 0.0;
            ValueType fc2;
            scai::IndexType NK = 256; 
            scai::lama::DenseVector<ValueType> freqVec; 
            scai::lama::DenseVector<ValueType> kVec;         
    };
}
