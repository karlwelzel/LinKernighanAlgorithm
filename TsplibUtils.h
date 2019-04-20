//
// Created by Karl Welzel on 19/04/2019.
//

#ifndef LINKERNINGHANALGORITHM_TSPLIBUTILS_H
#define LINKERNINGHANALGORITHM_TSPLIBUTILS_H

#include <string>
#include <vector>
#include "Tour.h"


std::string trim(const std::string &str, const std::string &whitespace = " \t");


class TsplibProblem {
private:
    const char delimiter = ':';

    std::string name;
    std::string type;
    dimension_t dimension = 0;
    std::string edgeWeightType;
    std::string edgeWeightFormat;
    std::string nodeCoordType = "TWOD_COORDS";

    // if EDGE_WEIGHT_TYPE is EXPLICIT
    // The matrix described in the TSPLIB file
    std::vector<std::vector<distance_t>> matrix;

    // if EDGE_WEIGHT_TYPE is *_2D
    std::vector<std::vector<double>> coordinates; // Every entry is a 2d coordinate

    std::string interpretKeyword(const std::string &keyword, const std::string &value);

public:
    TsplibProblem();

    std::string readFile(std::ifstream &inputFile);


    const std::string &getName() const;

    const std::string &getType() const;

    dimension_t getDimension() const;


    distance_t dist(vertex_t i, vertex_t j) const;

    distance_t length(const Tour &tour) const;
};


class TsplibTour : public Tour {
private:
    const char delimiter = ':';

    std::string name = "";
    std::string type = "";
    dimension_t dimension = 0;

    std::string interpretKeyword(const std::string &keyword, const std::string &value);

public:
    std::string readFile(std::ifstream &inputFile);

    const std::string &getName() const;

    const std::string &getType() const;
};

#endif //LINKERNINGHANALGORITHM_TSPLIBUTILS_H
