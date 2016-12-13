#include "Utilities.h"
#include <fstream>

std::string autoFileName(std::string prefix){

    std::time_t t = std::time(NULL);
    std::tm* timePtr = localtime(&t);
    char buffer[256];
    std::strftime(buffer, sizeof(buffer), "_%y-%m-%d_%H.%M.%S.txt", timePtr);

    return prefix + buffer;
}



void print_table(std::string prefix, const Eigen::MatrixXd & table) {
    auto fileName = autoFileName(prefix);
    std::fstream outFile;
    outFile.open(fileName, std::fstream::out);
    if(!outFile.good())
        throw std::runtime_error("Could not open file " + fileName + " for output. Aborting");


    for(int i = 0; i < table.cols(); i++){
        outFile << table.col(i).transpose() << std::endl;
    }
}
