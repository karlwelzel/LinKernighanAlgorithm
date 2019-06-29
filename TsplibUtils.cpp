//
// Created by Karl Welzel on 19.04.2019.
//

#include <cmath>
#include <cstddef>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "TsplibUtils.h"

// Checks if str has the UPPERCASE_WITH_UNDERSCORES format, that all keywords have
bool isKeyword(const std::string &str) {
    return str.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ_") == std::string::npos;
}

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

TsplibProblem::TsplibProblem(bool storeAllDistances) : storeAllDistances(storeAllDistances) {
}

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
            edgeWeightFormat != "UPPER_DIAG_ROW" and edgeWeightFormat != "UPPER_ROW") {
            return "EDGE_WEIGHT_FORMAT must be one of: FULL_MATRIX, LOWER_DIAG_ROW, UPPER_DIAG_ROW, UPPER_ROW";
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
    std::string::size_type delimiterIndex;
    std::vector<distance_t> numbers; // All numbers in the TSPLIB file as they appear in EDGE_WEIGHT_SECTION

    while (inputFile) {
        getline(inputFile, line);
        line = trim(line);
        if (line.empty()) {
            // ignore empty lines
        } else if ((delimiterIndex = line.find(DELIMITER)) != std::string::npos) {
            // The specification part with entries of the form <keyword> : <value>
            const std::string keyword = trim(line.substr(0, delimiterIndex));
            const std::string value = trim(line.substr(delimiterIndex + 1, std::string::npos));
            const std::string errorMessage = interpretKeyword(keyword, value);
            if (!errorMessage.empty()) {
                return errorMessage;
            }
        } else {
            // The data part or an invalid part of the file
            if (isKeyword(line)) {
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
    } else if (storeAllDistances) {
        matrix.assign(dimension, std::vector<distance_t>(dimension, 0));
        for (vertex_t i = 0; i < dimension; ++i) {
            for (vertex_t j = 0; j < dimension; ++j) {
                matrix[i][j] = trueDistance(i, j);
            }
        }
    }

    return "";
}

const std::string &TsplibProblem::getName() const {
    return name;
}

const std::string &TsplibProblem::getType() const {
    return type;
}

dimension_t TsplibProblem::getDimension() const {
    return dimension;
}

distance_t TsplibProblem::trueDistance(vertex_t i, vertex_t j) const {
    if (edgeWeightType == "EUC_2D") {
        double d = hypot(coordinates[i][0] - coordinates[j][0],
                         coordinates[i][1] - coordinates[j][1]);
        return static_cast<distance_t>(lround(d)); // lround(d) >= 0, so casting does not produce any errors
    } else if (edgeWeightType == "CEIL_2D") {
        double d = hypot(coordinates[i][0] - coordinates[j][0],
                         coordinates[i][1] - coordinates[j][1]);
        // ceil(d) returns a double and to prevent errors the result is rounded before casting to distance_t
        return static_cast<distance_t>(lround(ceil(d)));
    } else if (edgeWeightType == "EXPLICIT") {
        return matrix[i][j];
    } else {
        throw std::runtime_error("The EDGE_WEIGHT_TYPE '" + edgeWeightType + "' is not supported.");
    }
}

distance_t TsplibProblem::dist(const vertex_t i, const vertex_t j) const {
    if (storeAllDistances) {
        return matrix[i][j];
    } else {
        return trueDistance(i, j);
    }
}

distance_t TsplibProblem::length(const BaseTour &tour) const {
    distance_t sum = 0;
    vertex_t currentVertex = 0;
    vertex_t nextVertex = tour.successor(currentVertex);
    do {
        sum += dist(currentVertex, nextVertex);
        currentVertex = nextVertex;
        nextVertex = tour.successor(currentVertex);
    } while (currentVertex != 0);
    return sum;
}

signed_distance_t TsplibProblem::exchangeGain(std::vector<vertex_t> &alternatingWalk) const {
    signed_distance_t value = 0;
    for (size_t i = 0; i < alternatingWalk.size() - 1; ++i) {
        if (i % 2 == 0) {
            value += dist(alternatingWalk[i], alternatingWalk[i + 1]);
        } else {
            value -= dist(alternatingWalk[i], alternatingWalk[i + 1]);
        }
    }
    return value;
}


// ============================================== TsplibTour class =====================================================

TsplibTour::TsplibTour() = default;

TsplibTour::TsplibTour(std::string name, const Tour &tour) : Tour(tour), name(std::move(name)), type("TOUR"),
                                                             dimension(tour.getDimension()) {}

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

std::string TsplibTour::readFile(std::ifstream &inputFile) {
    std::string line;
    std::string lastDataKeyword;
    std::string::size_type delimiterIndex;
    std::vector<vertex_t> tourList;

    while (inputFile) {
        getline(inputFile, line);
        line = trim(line);
        if (line.empty()) {
            // ignore this line
        } else if ((delimiterIndex = line.find(DELIMITER)) != std::string::npos) {
            // The specification part with entries of the form <keyword> : <value>
            const std::string keyword = trim(line.substr(0, delimiterIndex));
            const std::string value = trim(line.substr(delimiterIndex + 1, std::string::npos));
            const std::string errorMessage = interpretKeyword(keyword, value);
            if (!errorMessage.empty()) {
                return errorMessage;
            }
        } else {
            // The data part or an invalid part of the file
            if (isKeyword(line)) {
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
    } else if (dimension != tourList.size()) {
        return "The dimension does not fit to the number of vertices";
    }
    setVertices(tourList);

    return "";
}

std::string TsplibTour::toTsplibTourFile() {
    std::string output;
    output += "NAME : " + name + "\n";
    output += "TYPE : " + type + "\n";
    output += "DIMENSION : " + std::to_string(dimension) + "\n";
    output += "TOUR_SECTION\n";
    vertex_t currentVertex = 0;
    do {
        output += std::to_string(currentVertex + 1) + "\n";
        currentVertex = successor(currentVertex);
    } while (currentVertex != 0);
    output += "-1\n";
    return output;
}

const std::string &TsplibTour::getName() const {
    return name;
}

const std::string &TsplibTour::getType() const {
    return type;
}


