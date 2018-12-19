#include <scai/lama.hpp>

#include "StepLengthSearch.hpp"
#include <gtest/gtest.h>

using namespace scai;
using namespace KITGPI;

typedef double ValueType;
 
bool verbose; // global variable definition
 
TEST(StepLengthSearchTest, TestParabolicFit)
{
lama::DenseVector<ValueType> xValues(3,1,scai::hmemo::Context::getContextPtr());
lama::DenseVector<ValueType> yValues(3,1,scai::hmemo::Context::getContextPtr());

xValues.setValue(0,-6);
xValues.setValue(1,-2);
xValues.setValue(2,2);
yValues.setValue(0,4);
yValues.setValue(1,2);
yValues.setValue(2,4);

ValueType solution=-2;

StepLengthSearch<ValueType> SLsearch;
ValueType res=SLsearch.parabolicFit(xValues,yValues);
EXPECT_EQ(solution,res);

lama::DenseVector<ValueType> testValues(4,1,scai::hmemo::Context::getContextPtr());
EXPECT_ANY_THROW(SLsearch.parabolicFit(testValues,yValues));
EXPECT_ANY_THROW(SLsearch.parabolicFit(xValues,testValues));
}


