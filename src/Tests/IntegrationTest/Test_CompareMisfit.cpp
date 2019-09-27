#include <iostream>
#include <scai/lama.hpp>
#include <string>
#include <vector>

#include "../Common/HostPrint.hpp"
#include <Configuration/Configuration.hpp>
#include <Configuration/ValueType.hpp>

using namespace scai;
using namespace KITGPI;

void skipHeaderLines(std::ifstream &ifs)
{
    char firstChar = ifs.peek();
    while (firstChar == '#') {
        ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        firstChar = ifs.peek();
    }
}

int main(int argc, char *argv[])
{

    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    // read configuration parameter from file
    Configuration::Configuration config(argv[1]);
    std::string filenameRef = config.get<std::string>("LogFilename");
    std::string dimensionRef = config.get<std::string>("dimension");
    std::string equationTypeRef = config.get<std::string>("equationType");

    // read reference misfit
    std::string filenameMisfit = "ci/ReferenceMisfits.txt";
    std::ifstream misfitFile;
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
    ValueType referenceMisfit(std::stof(referenceMisfitStr));

    //get smallest misfit value from log file
    std::ifstream logFile;
    logFile.open(filenameRef);
    skipHeaderLines(logFile);

    std::string bufferStr;
    std::string misfit;

    while (!logFile.eof() && !logFile.bad()) {
        logFile >> bufferStr >> bufferStr >> bufferStr >> bufferStr >> bufferStr >> bufferStr >> bufferStr >> bufferStr >> bufferStr >> bufferStr >> misfit;
    }

    ValueType finalMisfit = std::stof(misfit);
    logFile.close();

    if (finalMisfit > referenceMisfit) {
        std::cout << "\n\nTest Failed\n\n"
                  << std::endl
                  << "Inversion misfit of " << finalMisfit << " is greater than reference misfit of " << referenceMisfit << "\n\n"
                  << std::endl;
        return (1);
    } else {
        std::cout << "\n\n!!! Successful !!!\n\n"
                  << std::endl
                  << "Inversion misfit of " << finalMisfit << " is smaller than reference misfit of " << referenceMisfit << "\n\n"
                  << std::endl;
    }

    return 0;
}
