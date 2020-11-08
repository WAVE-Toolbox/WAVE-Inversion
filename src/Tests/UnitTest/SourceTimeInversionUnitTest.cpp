#include "../../SourceEstimation/SourceEstimation.hpp"
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <gtest/gtest.h>
#include <scai/lama.hpp>

using namespace scai;
using namespace KITGPI;
typedef double ValueType;

TEST(SourceTimeInversionTest, TestSourceEstimation)
{

    SourceEstimation<ValueType> sourceEst;

    dmemo::DistributionPtr dist(new dmemo::NoDistribution(10000));
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    Configuration::Configuration testConfig("../src/Tests/Testfiles/testSourceTimeInversion_config.txt");

    // Create an object of the mapping (3D-1D) class Coordinates

    Acquisition::Coordinates<ValueType> modelCoordinates(testConfig.get<IndexType>("NX"), testConfig.get<IndexType>("NY"), testConfig.get<IndexType>("NZ"), testConfig.get<ValueType>("DH"));

    Acquisition::Receivers<ValueType> receivers;
    receivers.init(testConfig, modelCoordinates, ctx, dist);
    lama::DenseMatrix<ValueType> &receiversData = receivers.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();
    receiversData.readFromFile("../src/Tests/Testfiles/testSourceTimeInversion_synth.shot_0.p.mtx");

    Acquisition::Receivers<ValueType> receiversTrue;
    receiversTrue.init(testConfig, modelCoordinates, ctx, dist);
    lama::DenseMatrix<ValueType> &receiversTrueData = receiversTrue.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();
    receiversTrueData.readFromFile("../src/Tests/Testfiles/testSourceTimeInversion_true.shot_0.p.mtx");

    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings;
    Acquisition::readAllSettings<ValueType>(sourceSettings, testConfig.get<std::string>("SourceFilename") + ".txt");
    Acquisition::Sources<ValueType> sources;
    sources.init(sourceSettings, testConfig, modelCoordinates, ctx, dist);

    lama::DenseMatrix<ValueType> &sourcesData = sources.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();
    sourcesData.readFromFile("../src/Tests/Testfiles/testSourceTimeInversion_sourceSignal.mtx");

    sourceEst.init(500, sources.get1DCoordinates().getDistributionPtr(), 1.0e-10);
    std::vector<scai::IndexType> filterHistoryCount(1, 0);
    sourceEst.estimateSourceSignal(receivers, receiversTrue, 0, 0, filterHistoryCount, useStreamConfig);
    sourceEst.applyFilter(sources, 0);

    lama::DenseMatrix<ValueType> sourceSignalInv = sources.getSeismogramHandler().getSeismogram(Acquisition::SeismogramType::P).getData();

    lama::DenseMatrix<ValueType> sourceSignalRef;
    sourceSignalRef.readFromFile("../src/Tests/Testfiles/testSourceTimeInversion_referenceSourceSignal.mtx");

    sourceSignalRef -= sourceSignalInv;
    EXPECT_LT(sourceSignalRef.l2Norm(), 0.4);
}
