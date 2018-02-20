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

    
    
    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    int width=40;
    int height=40;	

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

    
   //top-left
  for (IndexType x = NX/2-width/2; x <= NX/2; ++x) {
    for (IndexType y = NY/2-height/2; y <= NY/2; ++y) {
        vp(lama::Range(), y, x) *= 1.1;
    }
   }

    //top-right
  for (IndexType x = NX/2+1; x < NX/2+width/2; ++x) {
    for (IndexType y = NY/2-height/2; y <= NY/2; ++y) {
        vp(lama::Range(), y, x) *= 0.95;
    }
   }      

 //bottom left        
 for (IndexType x = NX/2-width/2; x <= NX/2; ++x) {
  for (IndexType y = NY/2+1; y < NY/2+width/2; ++y) {
      vp(lama::Range(), y, x) *= 0.9;
  }                                                               
 }

    //bottom right
  for (IndexType x = NX/2+1; x < NX/2+width/2; ++x) {
    for (IndexType y = NY/2+1; y < NY/2+width/2; ++y) {
        vp(lama::Range(), y, x) *= 1.05;
    }
   }      

    //write model to file specified in configuration
    std::string filename = config.get<std::string>("ModelFilename");
    vp.writeToFile(filename + ".vp.mtx");
    rho.writeToFile(filename + ".density.mtx");
    return 0;
}
