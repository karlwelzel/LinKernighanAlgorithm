//
// Created by Karl Welzel on 29/03/2019.
//

// This library contains code for handling TSP problems from the TSPLIB project

#include <utility>
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <regex>
#include "TsplibUtils.h"


// Trim whitespaces from the start and end of str
std::string trim(const std::string &str, const std::string &whitespace = " \t") {
    unsigned int start = str.find_first_not_of(whitespace);
    unsigned int end = str.find_last_not_of(whitespace);
    if (start == std::string::npos) {
        return ""; // No non-whitespace character, trimming leaves only the empty string
    } else {
        return str.substr(start, end - start + 1);
    }
}


// ============================================ TsplibProblem class ====================================================

std::string TsplibProblem::interpretKeyword(const std::string &keyword, const std::string &value) {
    if (keyword == "NAME") {
        name = value;
    } else if (keyword == "TYPE") {
        type = value;
        if (value != "TSP") {
            return "TYPE must be TSP";
        }

    } else if (keyword == "COMMENT") {
    } else if (keyword == "DIMENSION") {
        dimension = static_cast<unsigned int>(stoi(value));
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

std::string TsplibProblem::readFile(std::ifstream &inputFile) {
    std::string line;
    std::string lastDataKeyword;
    unsigned int delimiterIndex;
    std::vector<int> numbers; // All numbers in the TSPLIB file as they appear in EDGE_WEIGHT_SECTION

    while (inputFile) {
        getline(inputFile, line);
        line = trim(line);
        if (line.empty()) {
            // ignore this line
        } else if ((delimiterIndex = line.find(delimiter)) != std::string::npos) {
            // The specification part with entries of the form <keyword> : <value>
            std::string keyword = trim(line.substr(0, delimiterIndex));
            std::string value = trim(line.substr(delimiterIndex + 1, std::string::npos));
            std::string errorMessage = interpretKeyword(keyword, value);
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
                        return "NODE_COORD_SECTION encountered, but EDGE_WEIGHT_TYPE is not on of: EUC_2D, EUC_3D, "
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
                    unsigned int n = 0;
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
                    int n;
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
                return "Too few coordinates were specified or the coordinates have the wrong dimension (3D instead of 2D)";
            }
        }
    }

    // Fill the matrix of distances properly
    if (edgeWeightType == "EXPLICIT") {
        // Initialize the matrix with zeros
        matrix = std::vector<std::vector<int>>(dimension, std::vector<int>(dimension));
        try {
            unsigned int numbersIndex = 0;
            if (edgeWeightFormat == "FULL_MATRIX") {
                for (int i = 0; i < dimension; ++i) {
                    for (int j = 0; j < dimension; ++j) {
                        matrix[i][j] = numbers.at(numbersIndex);
                        numbersIndex++;
                    }
                }
            } else if (edgeWeightFormat == "LOWER_DIAG_ROW") {
                for (int i = 0; i < dimension; ++i) {
                    for (int j = 0; j <= i; ++j) {
                        matrix[i][j] = numbers.at(numbersIndex);
                        matrix[j][i] = numbers.at(numbersIndex);
                        numbersIndex++;
                    }
                }
            } else if (edgeWeightFormat == "UPPER_DIAG_ROW") {
                for (int i = 0; i < dimension; ++i) {
                    for (int j = i; j < dimension; ++j) {
                        matrix[i][j] = numbers.at(numbersIndex);
                        matrix[j][i] = numbers.at(numbersIndex);
                        numbersIndex++;
                    }
                }
            } else if (edgeWeightFormat == "UPPER_ROW") {
                // The diagonal is never touched, so it is filled with zeros from the initialization
                for (int i = 0; i < dimension; ++i) {
                    for (int j = i + 1; j < dimension; ++j) {
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

TsplibProblem::TsplibProblem() = default;

const std::string &TsplibProblem::getName() const {
    return name;
}

const std::string &TsplibProblem::getType() const {
    return type;
}

unsigned int TsplibProblem::getDimension() const {
    return dimension;
}

int TsplibProblem::dist(unsigned int i, unsigned int j) const {
    if (edgeWeightType == "EUC_2D") {
        double d = std::hypot(coordinates.at(i).at(0) - coordinates.at(j).at(0),
                              coordinates.at(i).at(1) - coordinates.at(j).at(1));
        return std::lround(d);
    } else if (edgeWeightType == "CEIL_2D") {
        double d = std::hypot(coordinates.at(i).at(0) - coordinates.at(j).at(0),
                              coordinates.at(i).at(1) - coordinates.at(j).at(1));
        return std::ceil(d);
    } else if (edgeWeightType == "EXPLICIT") {
        return matrix.at(i).at(j);
    } else {
        throw std::runtime_error("The EDGE_WEIGHT_TYPE '" + edgeWeightType + "' is not supported.");
    }
}


// =================================================== Tour class ======================================================

Tour::Tour() = default;

Tour::Tour(const std::list<unsigned int> &tour) : vertices(tour) {}

const std::list<unsigned int> &Tour::getVertices() const {
    return vertices;
}

const unsigned int Tour::length(TsplibProblem &tsplibProblem) {
    unsigned int sum = 0;
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        if (std::next(it) == vertices.end()) {
            sum += tsplibProblem.dist(*it, vertices.front());
        } else {
            sum += tsplibProblem.dist(*it, *std::next(it));
        }
    }
    return sum;
}

std::ostream &operator<<(std::ostream &out, const Tour &tour) {
    std::string output;
    for (unsigned int vertex : tour.getVertices()) {
        output += std::to_string(vertex + 1) + ", ";
    }
    out << output.substr(0, output.length() - 2) << std::endl;
    return out;
}


// ============================================== TsplibTour class =====================================================

const std::string &TsplibTour::getName() const {
    return name;
}

const std::string &TsplibTour::getType() const {
    return type;
}

std::string TsplibTour::interpretKeyword(const std::string &keyword, const std::string &value) {
    if (keyword == "NAME") {
        name = value;
    } else if (keyword == "TYPE") {
        type = value;
        if (value != "TOUR") {
            return "TYPE must be TOUR";
        }

    } else if (keyword == "COMMENT") {
    } else if (keyword == "DIMENSION") {
        dimension = static_cast<unsigned int>(stoi(value));
    } // every unknown keyword is ignored
    return "";
}

std::string TsplibTour::readFile(std::ifstream &inputFile) {
    std::string line;
    std::string lastDataKeyword;
    unsigned int delimiterIndex;

    while (inputFile) {
        getline(inputFile, line);
        line = trim(line);
        if (line.empty()) {
            // ignore this line
        } else if ((delimiterIndex = line.find(delimiter)) != std::string::npos) {
            // The specification part with entries of the form <keyword> : <value>
            std::string keyword = trim(line.substr(0, delimiterIndex));
            std::string value = trim(line.substr(delimiterIndex + 1, std::string::npos));
            std::string errorMessage = interpretKeyword(keyword, value);
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
                    int n;
                    while (stream >> n) {
                        if (n >= 1) {
                            vertices.push_back(n - 1);
                        }
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
    if (vertices.size() != dimension) {
        return "The dimension does not fit to the number of vertices";
    }

    return "";
}

