//
// Created by karl on 19.04.19.
//

#ifndef LINKERNINGHANALGORITHM_TSPLIBTOUR_H
#define LINKERNINGHANALGORITHM_TSPLIBTOUR_H

#include "Tour.h"

class TsplibTour : public Tour {
private:
    const char delimiter = ':';

    std::string name = "";
    std::string type = "";
    unsigned int dimension = 0;

    std::string interpretKeyword(const std::string &keyword, const std::string &value);

public:
    const std::string &getName() const;

    const std::string &getType() const;

    std::string readFile(std::ifstream &inputFile);
};

#endif //LINKERNINGHANALGORITHM_TSPLIBTOUR_H
