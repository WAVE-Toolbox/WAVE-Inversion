#include <scai/lama.hpp>
#include "../../SourceEstimation/SourceEstimation.hpp"
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <gtest/gtest.h>

using namespace scai;
using namespace KITGPI;
typedef double ValueType;

TEST(SourceTimeInversionTest, TestSourceEstimation)
{

    SourceEstimation<ValueType> sourceEst(0.01, 500);
    
    dmemo::DistributionPtr dist(new dmemo::NoDistribution(10000));
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();
    
    Configuration::Configuration testConfig("../src/Tests/Testfiles/testSourceTimeInversion_config.txt");
    
    Acquisition::Receivers<ValueType> receivers;
    receivers.init(testConfig, ctx, dist);
    lama::DenseMatrix<ValueType> &receiversData = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();
    receiversData.readFromFile("../src/Tests/Testfiles/testSourceTimeInversion_synth.shot_0.p.mtx");
    
    Acquisition::Receivers<ValueType> receiversTrue;
    receiversTrue.init(testConfig, ctx, dist);
    lama::DenseMatrix<ValueType> &receiversTrueData = receiversTrue.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();
    receiversTrueData.readFromFile("../src/Tests/Testfiles/testSourceTimeInversion_true.shot_0.p.mtx");
    
    Acquisition::Sources<ValueType> sources;
    sources.init(testConfig, ctx, dist);
    lama::DenseMatrix<ValueType> &sourcesData = sources.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();
    sourcesData.readFromFile("../src/Tests/Testfiles/testSourceTimeInversion_sourceSignal.mtx");
    
    sourceEst.estimateSourceSignal(receivers, receiversTrue, sources);
            
    lama::DenseMatrix<ValueType> sourceSignalInv = sources.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();
    
    lama::DenseMatrix<ValueType> sourceSignalRef;
    sourceSignalRef.readFromFile("../src/Tests/Testfiles/testSourceTimeInversion_referenceSourceSignal.mtx");
    
    sourceSignalRef -= sourceSignalInv;
    EXPECT_LT(sourceSignalRef.l2Norm(), 0.4);
    
}