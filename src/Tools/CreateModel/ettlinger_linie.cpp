#include <iostream>
#include <scai/common/Settings.hpp>
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/lama.hpp>
#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <string>
#include <vector>

#include <IO/IO.hpp>
#include <Configuration/Configuration.hpp>
#include "HostPrint.hpp"
#include <Configuration/ValueType.hpp>

using namespace KITGPI;
using namespace scai;
/*------------------------------
     Ettlinger Linie - Model Creation
     This is a simple model creation tool which
     reads in the configuration file the model dimension and name
-------------------------------*/
int main(int argc, char *argv[])
{
    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }
    
    // read configuration parameter from file
    KITGPI::Configuration::Configuration config(argv[1]);    
    
    // parameter gradient
    ValueType vp1_up=config.get<ValueType>("vp1_up");
    ValueType vp1_down=config.get<ValueType>("vp1_down");
    
    ValueType vs1_up=config.get<ValueType>("vs1_up");
    ValueType vs1_middle=config.get<ValueType>("vs1_middle");
    ValueType vs1_grad_change=config.get<ValueType>("depth_grad_change_vs");    
    ValueType vs1_down=config.get<ValueType>("vs1_down");
    
    // depth of second layer
    IndexType depth=config.get<ValueType>("depth_halfspace");
    //velocities of second layer
    ValueType vp2=config.get<ValueType>("vp2");
    ValueType vs2=config.get<ValueType>("vs2");
    
    ValueType rho1=config.get<ValueType>("rho1");
    ValueType rho2=config.get<ValueType>("rho2");
    
    // trench
    ValueType vs_trench=config.get<ValueType>("trench_vs");
    IndexType vs_centerx=config.get<ValueType>("trench_centerx");
    IndexType vs_width=config.get<ValueType>("trench_width");

    // estimate grid with parameters out of the configuration
    IndexType NX = config.get<IndexType>("NX");
    IndexType NY = config.get<IndexType>("NY");
    IndexType NZ = config.get<IndexType>("NZ");
    common::Grid3D grid(NY, NZ, NX);
    
    // construct model vectors and set velocities for first layer
    lama::GridVector<ValueType> vp(grid, vp1_up);
    lama::GridVector<ValueType> vs(grid, vs1_up);
    lama::GridVector<ValueType> rho(grid, rho1);
    
    // put gradients
    for (IndexType y = 0; y < depth; y++){
        ValueType temp = vp1_up + y*(vp1_down-vp1_up)/(depth-1);
        vp(y, lama::Range(), lama::Range()) = temp;
    }
    
    for (IndexType y = 0; y < vs1_grad_change; y++){
        ValueType temp = vs1_up + y*(vs1_middle-vs1_up)/(vs1_grad_change-1);
        vs(y, lama::Range(), lama::Range()) = temp;
    }
    
    for (IndexType y = vs1_grad_change; y < depth; y++){
        ValueType temp = vs1_middle + (y-vs1_grad_change)*(vs1_down-vs1_middle)/((depth-vs1_grad_change)-1);
        vs(y, lama::Range(), lama::Range()) = temp;
    }
    
    // set velocities for second layer
    for (IndexType y = depth; y < NY; ++y) {
        vp(y, lama::Range(), lama::Range()) = vp2;
        vs(y, lama::Range(), lama::Range()) = vs2;
        rho(y, lama::Range(), lama::Range()) = rho2;
    }
    
    // put trench
    for (IndexType x=0; x<NX; ++x) {
        for (IndexType y=0; y<NY; ++y) {
            if ((x>vs_centerx-2 && x+y<(vs_centerx-1)+vs_width/2+1) || (x<(vs_centerx-1)+1 && x-y>(vs_centerx-1)-vs_width/2-1)) {
                vs(y, lama::Range(), x) = vs_trench;
            }
        }
    }
      
    //write model to file specified in configuration
    std::string filename = config.get<std::string>("ModelFilename");
    
    IndexType fileFormat = config.get<IndexType>("FileFormat");
    
    //write model to disc
    KITGPI::IO::writeVector(vp, filename + ".vp", fileFormat);
    KITGPI::IO::writeVector(vs, filename + ".vs", fileFormat);
    KITGPI::IO::writeVector(rho, filename + ".density", fileFormat);

    return 0;
}
