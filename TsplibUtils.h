//
// Created by Karl Welzel on 19.04.2019.
//

#ifndef LINKERNIGHANALGORITHM_TSPLIBUTILS_H
#define LINKERNIGHANALGORITHM_TSPLIBUTILS_H

#include <fstream>
#include <string>
#include <vector>
#include "Tour.h"


std::string trim(const std::string &str, const std::string &whitespace = " \t");

// ============================================ TsplibProblem class ====================================================

// This class reads and stores problems from the TSPLIB library, while also checking for syntax errors, logical errors
// or unsupported keywords

class TsplibProblem {
private:
    // The delimiter between a keyword and its value, used in the TSPLIB format
    static const char DELIMITER = ':';

    std::string name;
    std::string type;
    dimension_t dimension = 0;
    std::string edgeWeightType;
    std::string edgeWeightFormat;
    std::string nodeCoordType = "TWOD_COORDS";

    // if EDGE_WEIGHT_TYPE is EXPLICIT or ALWAYS_STORE_IN_MATRIX is true
    // The matrix that contains all pairs of distances
    std::vector<std::vector<distance_t>> matrix;

    // if EDGE_WEIGHT_TYPE is *_2D
    std::vector<std::vector<double>> coordinates; // Every entry is a 2d coordinate

    // When the EDGE_WEIGHT_TYPE is not EXPLICIT, this parameter specifies whether all distances should be computed
    // and stored in a matrix or only the coordinates should be saved and the distances will be computed on every call
    // to TsplibProblem::dist. This is a trade-off between running time and memory usage
    bool storeAllDistances = true;

    // Interpret a single keyword value pair from the specification part of a TSPLIB file
    // Expects keyword and value to not have any superfluous whitespaces or line breaks
    // Returns an error message if an error occurred and an empty string otherwise
    std::string interpretKeyword(const std::string &keyword, const std::string &value);

    // Returns the distance of vertex i and vertex j and ignores ALWAYS_STORE_IN_MATRIX
    // Expects that i and j are in [0, dimension)
    distance_t trueDistance(vertex_t i, vertex_t j) const;

public:
    explicit TsplibProblem(bool storeAllDistances = true);

    // Interpret the file "inputFile" as a TSPLIB file and store the information given there
    // Returns an error message if an error occurred and an empty string otherwise
    std::string readFile(std::ifstream &inputFile);

    // Returns the name of the TSPLIB problem
    const std::string &getName() const;

    // Returns the type of the TSPLIB problem (always TSP)
    const std::string &getType() const;

    // Returns the number of vertices in the TSPLIB problem
    dimension_t getDimension() const;

    // Returns the distance of vertex i and vertex j
    // Expects that i and j are in [0, dimension)
    distance_t dist(vertex_t i, vertex_t j) const;

    // Returns the length of tour
    distance_t length(const BaseTour &tour) const;

    // Returns the gain of an alternating walk, i.e. the decrease in length of the tour if all the edges of the walk
    // on the tour are replaced by edges of the walk not on the tour
    signed_distance_t exchangeGain(const AlternatingWalk &alternatingWalk) const;
};


// ============================================== TsplibTour class =====================================================

// This class reads and stores tours from the TSPLIB library, while also checking for syntax errors, logical errors
// or unsupported keywords

class TsplibTour : public Tour {
private:
    const char DELIMITER = ':';

    std::string name = "";
    std::string type = "";
    dimension_t dimension = 0;

    // Interpret a single keyword value pair from the specification part of a TSPLIB tour file
    // Expects keyword and value to not have any superfluous whitespaces or line breaks
    // Returns an error message if an error occurred and an empty string otherwise
    std::string interpretKeyword(const std::string &keyword, const std::string &value);

public:
    TsplibTour();

    // Construct a TsplibTour from the tour and the name
    TsplibTour(std::string name, const Tour &tour);

    // Interpret the file "inputFile" as a TSPLIB tour file and store the information given there
    // Returns an error message if an error occurred and an empty string otherwise
    std::string readFile(std::ifstream &inputFile);

    // Converts the TsplibTour to the TSPLIB tour format
    std::string toTsplibTourFile();

    // Returns the name of the TSPLIB tour
    const std::string &getName() const;

    // Returns the type of the TSPLIB problem (always TOUR)
    const std::string &getType() const;
};

#endif //LINKERNIGHANALGORITHM_TSPLIBUTILS_H
