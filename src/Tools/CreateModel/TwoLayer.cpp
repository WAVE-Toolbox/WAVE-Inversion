#include <iostream>
#include <scai/common/Settings.hpp>
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/lama.hpp>
#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <string>
#include <vector>

#include "Configuration/Configuration.hpp"
#include "HostPrint.hpp"

using namespace scai;

int main(int argc, char *argv[])
{
    typedef double ValueType;

    // parameter contrast for the second layer
    ValueType contrast = 1.3;
    
    
    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    // read configuration parameter from file
    KITGPI::Configuration::Configuration config(argv[1]);

    // estimate grid with parameters out of the configuration
    int NX = config.get<IndexType>("NX");
    int NY = config.get<IndexType>("NY");
    int NZ = config.get<IndexType>("NZ");
    common::Grid3D grid(NZ, NY, NX);

    // construct model vectors
    lama::GridVector<ValueType> vp(grid);
    lama::GridVector<ValueType> vs(grid);
    lama::GridVector<ValueType> rho(grid);

    //set models to values specified in configuration
    vp = config.get<ValueType>("velocityP");;
    vs = config.get<ValueType>("velocityS");
    rho = config.get<ValueType>("rho");

    
    //set velocities for second layer
    for (IndexType y = NY / 2; y < NY; ++y) {

        vp(lama::Range(), y, lama::Range()) *= contrast;
        vs(lama::Range(), y, lama::Range()) *= contrast;
        rho(lama::Range(), y, lama::Range()) *= contrast;
    }
    
    //write model to file specified in configuration
    std::string filename = config.get<std::string>("ModelFilename");
    vp.writeToFile(filename + ".vp.mtx");
    vs.writeToFile(filename + ".vs.mtx");
    rho.writeToFile(filename + ".density.mtx");
    return 0;
}
