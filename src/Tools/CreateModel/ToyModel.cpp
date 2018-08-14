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

    ValueType vp1, vp2, vp3, vp4, vs1, vs2, vs3, vs4, density1, density2, density3, density4;
    int depth=1.0;
    
    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }
    // read configuration parameter from file
    KITGPI::Configuration::Configuration config(argv[1]);
    
    //background value
    ValueType vp0=config.get<ValueType>("velocityP");
    ValueType vs0=config.get<ValueType>("velocityS");
    ValueType density0=config.get<ValueType>("rho");
    
     // box
    int width=config.get<IndexType>("boxWidth");
    int height=config.get<IndexType>("boxHeight");
    
    if (config.get<std::string>("dimension")=="3D")
    depth=config.get<IndexType>("boxDepth");

//    permutations in box
    try {vp1=config.get<ValueType>("vp1");} catch (...) { vp1=vp0;};
    try {vp2=config.get<ValueType>("vp2");} catch (...) { vp2=vp0;};
    try {vp3=config.get<ValueType>("vp3");} catch (...) { vp3=vp0;};
    try {vp4=config.get<ValueType>("vp4");} catch (...) { vp4=vp0;};

    try {vs1=config.get<ValueType>("vs1");} catch (...) { vs1=vs0;};
    try {vs2=config.get<ValueType>("vs2");} catch (...) { vs2=vs0;};
    try {vs3=config.get<ValueType>("vs3");} catch (...) { vs3=vs0;};
    try {vs4=config.get<ValueType>("vs4");} catch (...) { vs4=vs0;};
    
    try {density1=config.get<ValueType>("density1");} catch (...) { density1=density0;};
    try {density2=config.get<ValueType>("density2");} catch (...) { density2=density0;};
    try {density3=config.get<ValueType>("density3");} catch (...) { density3=density0;};
    try {density4=config.get<ValueType>("density4");} catch (...) { density4=density0;};

    // estimate grid with parameters out of the configuration
    int NX = config.get<IndexType>("NX");
    int NY = config.get<IndexType>("NY");
    int NZ = config.get<IndexType>("NZ");
    common::Grid3D grid(NZ, NY, NX);

    // construct model vectors
    lama::GridVector<ValueType> vp(grid);
    lama::GridVector<ValueType> vs(grid);
    lama::GridVector<ValueType> rho(grid);
    lama::GridVector<ValueType> tauP(grid);
    lama::GridVector<ValueType> tauS(grid);
    

    //set background value
    vp = vp0;
    vs = vs0;
    rho = density0;
    tauP = config.get<ValueType>("tauP");
    tauS = config.get<ValueType>("tauS");
    
    
   //top-left
  for (IndexType x = NX/2-width/2; x < NX/2; ++x) {
    for (IndexType y = NY/2-height/2; y < NY/2; ++y) {
	for (IndexType z = NZ/2-depth/2; z <= NZ/2+depth/2; ++z) {  
	  vp(z, y, x) = vp1; 
	  vs(z, y, x) = vs1; 
          rho(z, y, x) = density1;
	}
    }
   }

    //top-right
  for (IndexType x = NX/2; x < NX/2+width/2; ++x) {
    for (IndexType y = NY/2-height/2; y < NY/2; ++y) {    
	for (IndexType z = NZ/2-depth/2; z <= NZ/2+depth/2; ++z) {      
	  vp(z, y, x) = vp2; 
	  vs(z, y, x) = vs2; 
          rho(z, y, x) = density2;
	}
    }
   }      

 //bottom left        
 for (IndexType x = NX/2-width/2; x < NX/2; ++x) {
  for (IndexType y = NY/2; y < NY/2+width/2; ++y) {
      for (IndexType z = NZ/2-depth/2; z <= NZ/2+depth/2; ++z) {  
	  vp(z, y, x) = vp3; 
	  vs(z, y, x) = vs3; 
          rho(z, y, x) = density3;
      }
  }                                                               
 }

    //bottom right
  for (IndexType x = NX/2; x < NX/2+width/2; ++x) {
    for (IndexType y = NY/2; y < NY/2+width/2; ++y) {
	for (IndexType z = NZ/2-depth/2; z <= NZ/2+depth/2; ++z) {  
	  vp(z, y, x) = vp4; 
	  vs(z, y, x) = vs4; 
          rho(z, y, x) = density4;
	}
    }
   }      

   std::string type = config.get<std::string>("equationType");
   
    //write model to file specified in configuration
    std::string filename = config.get<std::string>("ModelFilename");
    
    rho.writeToFile(filename + ".density.mtx");

    if (type.compare("sh") != 0) {
        vp.writeToFile(filename + ".vp.mtx");
    }

    if (type.compare("acoustic") != 0) {
        vs.writeToFile(filename + ".vs.mtx");
    }
    if (type.compare("visco") == 0) {
        tauP.writeToFile(filename + ".tauP.mtx");
        tauS.writeToFile(filename + ".tauS.mtx");
    }
    return 0;
}
