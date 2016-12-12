#include "Utilities.h"

std::string autoFileName(){

    std::time_t t = std::time(NULL);
    std::tm* timePtr = localtime(&t);
    char buffer[256];
    std::strftime(buffer, sizeof(buffer), "output_%Y-%m-%d_%H.%M.%S.txt", timePtr);

    return buffer;
}
