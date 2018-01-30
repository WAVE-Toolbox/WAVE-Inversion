#include <iostream>
#include <vector>
#include <string>
#include <scai/lama.hpp>

#include "../Common/HostPrint.hpp"
#include <Configuration/Configuration.hpp>


void skipHeaderLines(std::ifstream &ifs) {
  char firstChar = ifs.peek();
  while (firstChar == '#') {
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    firstChar = ifs.peek();
  }
}

int main(int argc, char *argv[])
{
    typedef double ValueType;
    
    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    // read configuration parameter from file
    KITGPI::Configuration::Configuration config(argv[1]);
    std::string filenameRef = config.get<std::string>("LogFilename");
    std::string dimensionRef = config.get<std::string>("dimension");
    std::string equationTypeRef = config.get<std::string>("equationType");
    
    // read reference misfit
    std::string filenameMisfit = "ci/ReferenceMisfits.txt";
    std:: ifstream misfitFile;
    misfitFile.open(filenameMisfit);
    skipHeaderLines(misfitFile);
    
    bool flagCorrectLine = false;
    std::string dimension;
    std::string equationType;
    std::string referenceMisfitStr;

    while (!flagCorrectLine) {
      if (misfitFile.eof()) {
	std::cout << "\n\nNo reference misfit found for given dimension and equationType!\n\n"
                  << std::endl;
        return (2);
      }
      misfitFile >> dimension >> equationType >> referenceMisfitStr;
      if (dimension == dimensionRef && equationType == equationTypeRef) {
	flagCorrectLine = true;
      }
    }
    misfitFile.close();
    scai::lama::Scalar referenceMisfit(std::stof(referenceMisfitStr));
    
    //get smallest misfit value from log file
    std::ifstream logFile;    
    logFile.open(filenameRef);
    skipHeaderLines(logFile);
    
    std::string bufferStr;
    scai::lama::DenseVector<ValueType> misfits;
    misfits.allocate(3);
    std::string misfitfg, misfitsg, misfittg;
    scai::lama::Scalar iterationMin, misfitMin(referenceMisfit+1);
 
    while (!logFile.eof() && !logFile.bad()) {
      logFile >> bufferStr >> bufferStr >> bufferStr >> bufferStr >> bufferStr >> bufferStr >> misfitfg >> misfitsg >> misfittg >> bufferStr;     
      misfits.setValue(0,std::stof(misfitfg.substr(7,12)));
      misfits.setValue(1,std::stof(misfitsg.substr(7,12)));
      misfits.setValue(2,std::stof(misfittg.substr(7,12)));
      iterationMin = misfits.min();
      if (iterationMin < misfitMin) {
	misfitMin = iterationMin;
      }
    }
    
    logFile.close();

    if (misfitMin > referenceMisfit) {
        std::cout << "\n\nTest Failed\n\n"
                  << std::endl
                  << "Inversion misfit of " << misfitMin << " is smaller than reference misfit of " << referenceMisfit
                  << std::endl;
        return (1);
    } else {
        std::cout << "\n\n!!! Successful !!!\n\n"
                  << std::endl
                  << "Inversion misfit of " << misfitMin << " is smaller than reference misfit of " << referenceMisfit
                  << std::endl;
    }

    return 0;
}
