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

    ValueType vpTopLeftScale=1.0, vpTopRightScale=1.0, vpBottomLeftScale=1.0, vpBottomRightScale=1.0, vsTopLeftScale=1.0, vsTopRightScale=1.0, vsBottomLeftScale=1.0, vsBottomRightScale=1.0, rhoTopLeftScale=1.0, rhoTopRightScale=1.0, rhoBottomLeftScale=1.0, rhoBottomRightScale=1.0;
    
    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    // box
    int width=20;
    int height=20;	


//    permutations in box
    
    vpTopLeftScale=0.8;
    vpTopRightScale=1.0;
    vpBottomLeftScale=0.9;
    vpBottomRightScale=0.9;
//     
//     vsTopLeftScale=0.8;
//     vsTopRightScale=0.8;
//     vsBottomLeftScale=1.6;
//     vsBottomRightScale=1.6;
    
     rhoTopLeftScale=1.0;
     rhoTopRightScale=1.2;
     rhoBottomLeftScale=0.9;
     rhoBottomRightScale=0.9;

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
  for (IndexType x = NX/2-width/2; x < NX/2; ++x) {
    for (IndexType y = NY/2-height/2; y < NY/2; ++y) {
	  vp(lama::Range(), y, x) *= vpTopLeftScale; 
	  vs(lama::Range(), y, x) *= vsTopLeftScale; 
          rho(lama::Range(), y, x) *= rhoTopLeftScale;
    }
   }

    //top-right
  for (IndexType x = NX/2; x < NX/2+width/2; ++x) {
    for (IndexType y = NY/2-height/2; y < NY/2; ++y) {    
	  vp(lama::Range(), y, x) *= vpTopRightScale; 
	  vs(lama::Range(), y, x) *= vsTopRightScale; 
          rho(lama::Range(), y, x) *= rhoTopRightScale;
    }
   }      

 //bottom left        
 for (IndexType x = NX/2-width/2; x < NX/2; ++x) {
  for (IndexType y = NY/2; y < NY/2+width/2; ++y) {
	  vp(lama::Range(), y, x) *= vpBottomLeftScale; 
	  vs(lama::Range(), y, x) *= vsBottomLeftScale; 
          rho(lama::Range(), y, x) *= rhoBottomLeftScale;
  }                                                               
 }

    //bottom right
  for (IndexType x = NX/2; x < NX/2+width/2; ++x) {
    for (IndexType y = NY/2; y < NY/2+width/2; ++y) {
	  vp(lama::Range(), y, x) *= vpBottomRightScale; 
	  vs(lama::Range(), y, x) *= vsBottomRightScale; 
          rho(lama::Range(), y, x) *= rhoBottomRightScale;
    }
   }      

    //write model to file specified in configuration
    std::string filename = config.get<std::string>("ModelFilename");
    vp.writeToFile(filename + ".vp.mtx");
    vs.writeToFile(filename + ".vs.mtx");
    rho.writeToFile(filename + ".density.mtx");
    return 0;
}
