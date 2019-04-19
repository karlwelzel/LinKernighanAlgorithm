//
// Created by Karl Welzel on 19/04/2019.
//

#ifndef LINKERNINGHANALGORITHM_TSPLIBPROBLEM_H
#define LINKERNINGHANALGORITHM_TSPLIBPROBLEM_H

#include <string>
#include <vector>

std::string trim(const std::string &str, const std::string &whitespace = " \t");

class TsplibProblem {
private:
    const char delimiter = ':';

    std::string name;
    std::string type;
    unsigned int dimension = 0;
    std::string edgeWeightType;
    std::string edgeWeightFormat;
    std::string nodeCoordType = "TWOD_COORDS";

    // if EDGE_WEIGHT_TYPE is EXPLICIT
    std::vector<std::vector<int>, std::allocator<std::vector<int>>> matrix; // The matrix described in the TSPLIB file

    // if EDGE_WEIGHT_TYPE is *_2D
    std::vector<std::vector<double>> coordinates; // Every entry is a 2d coordinate

    std::string interpretKeyword(const std::string &keyword, const std::string &value);

public:
    TsplibProblem();

    std::string readFile(std::ifstream &inputFile);


    const std::string &getName() const;

    const std::string &getType() const;

    unsigned int getDimension() const;


    int dist(unsigned int i, unsigned int j) const;
};

#endif //LINKERNINGHANALGORITHM_TSPLIBPROBLEM_H
