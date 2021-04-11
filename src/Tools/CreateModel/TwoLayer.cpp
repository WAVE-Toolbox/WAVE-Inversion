#include <iostream>
#include <scai/common/Settings.hpp>
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/lama.hpp>
#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <string>
#include <vector>

#include <Configuration/Configuration.hpp>
#include <Configuration/ValueType.hpp>
#include "HostPrint.hpp"
#include <IO/IO.hpp>

using namespace scai;
using namespace KITGPI;

int main(int argc, char *argv[])
{

    // parameter contrast for the second layer
    ValueType contrast = 1.3;
    
    
    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    // read configuration parameter from file
    Configuration::Configuration config(argv[1]);

    // estimate grid with parameters out of the configuration
    int NX = config.get<IndexType>("NX");
    int NY = config.get<IndexType>("NY");
    int NZ = config.get<IndexType>("NZ");
    common::Grid3D grid(NY, NZ, NX);

    // construct model vectors and set value from configuration file
    lama::GridVector<ValueType> vp(grid,config.get<ValueType>("velocityP"));
    lama::GridVector<ValueType> vs(grid,config.get<ValueType>("velocityS"));
    lama::GridVector<ValueType> rho(grid,config.get<ValueType>("rho"));
    lama::GridVector<ValueType> tauP(grid,config.get<ValueType>("tauP"));
    lama::GridVector<ValueType> tauS(grid,config.get<ValueType>("tauS"));
    
    //set velocities for second layer
    for (IndexType y = NY / 2; y < NY; ++y) {

        vp( y, lama::Range(),lama::Range()) *= contrast;
        vs( y,lama::Range(), lama::Range()) *= contrast;
        rho(y,lama::Range(), lama::Range()) *= contrast;
    }
    
    
    std::string type = config.get<std::string>("equationType");
    std::transform(type.begin(), type.end(), type.begin(), ::tolower); 
    //write model to file specified in configuration
    std::string filename = config.get<std::string>("ModelFilename");
    
    IndexType fileFormat = config.get<IndexType>("FileFormat");

    //write model to disc

    KITGPI::IO::writeVector(rho, filename + ".density", fileFormat);

    if (type.compare("sh") != 0) {
        KITGPI::IO::writeVector(vp, filename + ".vp", fileFormat);
    }

    if (type.compare("acoustic") != 0) {
        KITGPI::IO::writeVector(vs, filename + ".vs", fileFormat);
    }

    if (type.compare("viscoelastic") == 0) {
        KITGPI::IO::writeVector(tauP, filename + ".tauP", fileFormat);
        KITGPI::IO::writeVector(tauS, filename + ".tauS", fileFormat);
    }
    return 0;
}
