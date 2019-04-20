//
// Created by Karl Welzel on 19/04/2019.
//

#include <list>
#include <cmath>
#include <regex>
#include <string>
#include <vector>
#include <fstream>
#include "TsplibUtils.h"


// Trim whitespaces from the start and end of str
std::string trim(const std::string &str, const std::string &whitespace) {
    const std::string::size_type start = str.find_first_not_of(whitespace);
    const std::string::size_type end = str.find_last_not_of(whitespace);
    if (start == std::string::npos) {
        return ""; // No non-whitespace character, trimming leaves only the empty string
    } else {
        return str.substr(start, end - start + 1);
    }
}


// ============================================ TsplibProblem class ====================================================

// This class reads and stores problems from the TSPLIB library, while also checking for syntax errors, logical errors
// or unsupported keywords


// Interpret a single keyword value pair from the specification part of a TSPLIB file
// Expects keyword and value to not have any superfluous whitespaces or line breaks
// Returns an error message if an error occurred and an empty string otherwise
std::string TsplibProblem::interpretKeyword(const std::string &keyword, const std::string &value) {
    if (keyword == "NAME") {
        name = value;
    } else if (keyword == "TYPE") {
        type = value;
        if (value != "TSP") {
            return "TYPE must be TSP";
        }
    } else if (keyword == "COMMENT") {
        // ignore any comments
    } else if (keyword == "DIMENSION") {
        dimension = static_cast<dimension_t>(stoi(value));
    } else if (keyword == "EDGE_WEIGHT_TYPE") {
        edgeWeightType = value;
        if (edgeWeightType != "EUC_2D" and edgeWeightType != "CEIL_2D" and edgeWeightType != "EXPLICIT") {
            return "EDGE_WEIGHT_TYPE must be one of: EUC_2D, CEIL_2D, EXPLICIT";
        }
    } else if (keyword == "EDGE_WEIGHT_FORMAT") {
        edgeWeightFormat = value;
        if (edgeWeightFormat != "FULL_MATRIX" and edgeWeightFormat != "LOWER_DIAG_ROW" and
            edgeWeightFormat != "UPPER_DIAG_ROW" and edgeWeightFormat != "UPPER_DIAG") {
            return "EDGE_WEIGHT_FORMAT must be one of: FULL_MATRIX, LOWER_DIAG_ROW, UPPER_DIAG_ROW, UPPER_DIAG";
        }
    } else if (keyword == "NODE_COORD_TYPE") {
        nodeCoordType = value;
        if (nodeCoordType != "TWOD_COORDS") {
            return "NODE_COORD_TYPE must be TWOD_COORDS";
        }
    } // every unknown keyword is ignored
    return "";
}

// Interpret the file "inputFile" as a TSPLIB file and store the information given there
// Returns an error message if an error occurred and an empty string otherwise
std::string TsplibProblem::readFile(std::ifstream &inputFile) {
    std::string line;
    std::string lastDataKeyword;
    std::string::size_type delimiterIndex;
    std::vector<distance_t> numbers; // All numbers in the TSPLIB file as they appear in EDGE_WEIGHT_SECTION

    while (inputFile) {
        getline(inputFile, line);
        line = trim(line);
        if (line.empty()) {
            // ignore empty lines
        } else if ((delimiterIndex = line.find(delimiter)) != std::string::npos) {
            // The specification part with entries of the form <keyword> : <value>
            const std::string keyword = trim(line.substr(0, delimiterIndex));
            const std::string value = trim(line.substr(delimiterIndex + 1, std::string::npos));
            const std::string errorMessage = interpretKeyword(keyword, value);
            if (!errorMessage.empty()) {
                return errorMessage;
            }
        } else {
            // The data part or an invalid part of the file
            if (std::regex_match(line, std::regex("[A-Z_]+"))) {
                // line is some kind of keyword, because all keywords have the UPPERCASE_WITH_UNDERSCORES format.
                lastDataKeyword = line;
                if (line == "EOF") {
                    break;
                } else if (line == "NODE_COORD_SECTION") {
                    if (edgeWeightType != "EUC_2D" and edgeWeightType != "EUC_3D" and edgeWeightType != "MAX_2D" and
                        edgeWeightType != "MAX_3D" and edgeWeightType != "MAN_2D" and edgeWeightType != "MAN_3D" and
                        edgeWeightType != "CEIL_2D" and edgeWeightType != "GEO") {
                        return "NODE_COORD_SECTION encountered, but EDGE_WEIGHT_TYPE is not one of: EUC_2D, EUC_3D, "
                               "MAX_2D, MAX_3D, MAN_2D, MAN_3D, CEIL_2D, GEO";
                    }
                    // Initialize the coordinates vector with empty vectors until coordinates.size() == dimension
                    while (coordinates.size() < dimension) {
                        coordinates.emplace_back();
                    }
                } else if (line == "EDGE_WEIGHT_SECTION") {
                    if (edgeWeightType != "EXPLICIT") {
                        return "EDGE_WEIGHT_SECTION encountered, but EDGE_WEIGHT_TYPE is not EXPLICIT";
                    }
                }
            } else {
                if (lastDataKeyword == "NODE_COORD_SECTION") {
                    // Already checked: nodeCoordType == "TWOD_COORDS"
                    // Every line has the format <integer> <real> <real>
                    std::stringstream stream(line);
                    vertex_t n = 0;
                    double x = 0, y = 0;
                    stream >> n >> x >> y;
                    try {
                        coordinates.at(n - 1) = std::vector<double>{x, y};
                    } catch (std::out_of_range &error) {
                        return "One of the coordinates has a number outside of the range [1, DIMENSION]";
                    }
                } else if (lastDataKeyword == "EDGE_WEIGHT_SECTION") {
                    // Every line is a sequence of integers separated by whitespaces
                    std::stringstream stream(line);
                    distance_t n;
                    while (stream >> n) {
                        numbers.push_back(n);
                    }
                } else if (!lastDataKeyword.empty()) {
                    // lastDataKeyword is unknown, so this block of data will be skipped
                } else {
                    return "Encountered data when expecting a keyword";
                }
            }
        }
    }

    // Error checking
    if (dimension == 0) {
        return "The dimension cannot be 0";
    }
    if (edgeWeightType.empty()) {
        return "EDGE_WEIGHT_TYPE must be specified";
    }

    if (edgeWeightType == "EXPLICIT" and edgeWeightFormat.empty()) {
        return "When EDGE_WEIGHT_TYPE is EXPLICIT, EDGE_WEIGHT_FORMAT must be specified";
    }

    if (edgeWeightType == "EUC_2D" or edgeWeightType == "MAX_2D" or edgeWeightType == "MAN_2D"
        or edgeWeightType == "CEIL_2D") {
        for (const std::vector<double> &coord : coordinates) {
            if (coord.size() != 2) {
                return "Too few coordinates were specified or the coordinates have the wrong dimension "
                       "(3D instead of 2D)";
            }
        }
    }

    // Fill the matrix of distances properly
    if (edgeWeightType == "EXPLICIT") {
        // Initialize the matrix with zeros
        matrix.assign(dimension, std::vector<distance_t>(dimension, 0));
        try {
            // TODO: Test all these possible edge weight formats
            size_t numbersIndex = 0;
            if (edgeWeightFormat == "FULL_MATRIX") {
                for (vertex_t i = 0; i < dimension; ++i) {
                    for (vertex_t j = 0; j < dimension; ++j) {
                        matrix[i][j] = numbers.at(numbersIndex);
                        numbersIndex++;
                    }
                }
            } else if (edgeWeightFormat == "LOWER_DIAG_ROW") {
                for (vertex_t i = 0; i < dimension; ++i) {
                    for (vertex_t j = 0; j <= i; ++j) {
                        matrix[i][j] = numbers.at(numbersIndex);
                        matrix[j][i] = numbers.at(numbersIndex);
                        numbersIndex++;
                    }
                }
            } else if (edgeWeightFormat == "UPPER_DIAG_ROW") {
                for (vertex_t i = 0; i < dimension; ++i) {
                    for (vertex_t j = i; j < dimension; ++j) {
                        matrix[i][j] = numbers.at(numbersIndex);
                        matrix[j][i] = numbers.at(numbersIndex);
                        numbersIndex++;
                    }
                }
            } else if (edgeWeightFormat == "UPPER_ROW") {
                // The diagonal is never touched, so it is filled with zeros from the initialization
                for (vertex_t i = 0; i < dimension; ++i) {
                    for (vertex_t j = i + 1; j < dimension; ++j) {
                        matrix[i][j] = numbers.at(numbersIndex);
                        matrix[j][i] = numbers.at(numbersIndex);
                        numbersIndex++;
                    }
                }
            }
            if (numbersIndex < numbers.size()) {
                return "Too many numbers were specified under EDGE_WEIGHT_SECTION";
            }
        } catch (std::out_of_range &error) {
            return "Too few numbers were specified under EDGE_WEIGHT_SECTION";
        }
    }

    return "";
}

// Returns the name of the TSPLIB problem
const std::string &TsplibProblem::getName() const {
    return name;
}

// Returns the type of the TSPLIB problem (always TSP)
const std::string &TsplibProblem::getType() const {
    return type;
}

// Returns the number of vertices in the TSPLIB problem
dimension_t TsplibProblem::getDimension() const {
    return dimension;
}

// Compute the distance of i and j
distance_t TsplibProblem::dist(const vertex_t i, const vertex_t j) const {
    if (edgeWeightType == "EUC_2D") {
        double d = hypot(coordinates.at(i).at(0) - coordinates.at(j).at(0),
                         coordinates.at(i).at(1) - coordinates.at(j).at(1));
        return static_cast<distance_t>(lround(d)); // lround(d) >= 0, so casting does not produce any
    } else if (edgeWeightType == "CEIL_2D") {
        double d = hypot(coordinates.at(i).at(0) - coordinates.at(j).at(0),
                         coordinates.at(i).at(1) - coordinates.at(j).at(1));
        // ceil(d) returns a double and to prevent errors the result is rounded before casting to distance_t
        return static_cast<distance_t>(lround(ceil(d)));
    } else if (edgeWeightType == "EXPLICIT") {
        return matrix.at(i).at(j);
    } else {
        throw std::runtime_error("The EDGE_WEIGHT_TYPE '" + edgeWeightType + "' is not supported.");
    }
}

distance_t TsplibProblem::length(const Tour &tour) const {
    distance_t sum = 0;
    TourWalker tourWalker(tour, 0);
    vertex_t currentVertex = tourWalker.getNextVertex();
    vertex_t nextVertex = tourWalker.getNextVertex();
    do {
        sum += dist(currentVertex, nextVertex);
        currentVertex = nextVertex;
        nextVertex = tourWalker.getNextVertex();
    } while (currentVertex != 0);
    return sum;
}


// ============================================== TsplibTour class =====================================================

// This class reads and stores tours from the TSPLIB library, while also checking for syntax errors, logical errors
// or unsupported keywords


TsplibProblem::TsplibProblem() = default;

// Interpret a single keyword value pair from the specification part of a TSPLIB tour file
// Expects keyword and value to not have any superfluous whitespaces or line breaks
// Returns an error message if an error occurred and an empty string otherwise
std::string TsplibTour::interpretKeyword(const std::string &keyword, const std::string &value) {
    if (keyword == "NAME") {
        name = value;
    } else if (keyword == "TYPE") {
        type = value;
        if (value != "TOUR") {
            return "TYPE must be TOUR";
        }

    } else if (keyword == "COMMENT") {
        // ignore any comments
    } else if (keyword == "DIMENSION") {
        dimension = static_cast<dimension_t>(stoi(value));
    } // every unknown keyword is ignored
    return "";
}

// Interpret the file "inputFile" as a TSPLIB tour file and store the information given there
// Returns an error message if an error occurred and an empty string otherwise
std::string TsplibTour::readFile(std::ifstream &inputFile) {
    std::string line;
    std::string lastDataKeyword;
    std::string::size_type delimiterIndex;
    std::list<vertex_t> tourList;

    while (inputFile) {
        getline(inputFile, line);
        line = trim(line);
        if (line.empty()) {
            // ignore this line
        } else if ((delimiterIndex = line.find(delimiter)) != std::string::npos) {
            // The specification part with entries of the form <keyword> : <value>
            const std::string keyword = trim(line.substr(0, delimiterIndex));
            const std::string value = trim(line.substr(delimiterIndex + 1, std::string::npos));
            const std::string errorMessage = interpretKeyword(keyword, value);
            if (!errorMessage.empty()) {
                return errorMessage;
            }
        } else {
            // The data part or an invalid part of the file
            if (std::regex_match(line, std::regex("[A-Z_]+"))) {
                // line is some kind of keyword, because all keywords have the UPPERCASE_WITH_UNDERSCORES format.
                lastDataKeyword = line;
                if (line == "EOF") {
                    break;
                }
            } else {
                if (lastDataKeyword == "TOUR_SECTION") {
                    // Every line is one single integer
                    std::stringstream stream(line);
                    vertex_t number;
                    if (line == "-1") {
                        lastDataKeyword.clear(); // The TOUR_SECTION is over
                    } else if (stream >> number and number >= 1) {
                        tourList.push_back(number - 1);
                    } else {
                        return "Expected a vertex or -1 in TOUR_SECTION, but got '" + line + "'";
                    }
                } else if (!lastDataKeyword.empty()) {
                    // lastDataKeyword is unknown, so this block of data will be skipped
                } else {
                    return "Encountered data when expecting a keyword";
                }
            }
        }
    }

    // Error checking
    if (dimension == 0) {
        return "The dimension cannot be 0";
    }
    setVertices(tourList);
    if (neighbors.size() != dimension) {
        return "The dimension does not fit to the number of vertices";
    }

    return "";
}

// Returns the name of the TSPLIB tour
const std::string &TsplibTour::getName() const {
    return name;
}

// Returns the type of the TSPLIB problem (always TOUR)
const std::string &TsplibTour::getType() const {
    return type;
}

