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
    std::string::size_type start = str.find_first_not_of(whitespace);
    std::string::size_type end = str.find_last_not_of(whitespace);
    if (start == std::string::npos) {
        return ""; // No non-whitespace character, trimming leaves only the empty string
    } else {
        return str.substr(start, end - start + 1);
    }
}


// ================================================ VertexList class ===================================================

// This class maintains a map that stores a maximum of two neighbors of each vertex. There is no inherent direction
// encoded in these neighbors. This gives additional flexibility when compared to std::list while still having similar
// run times. It is possible to create tours that do not contain every vertex or multiple separate paths or even have a
// vertex x with neighbor y while y does not have x as its neighbor, so extra caution is necessary when manipulating
// these objects.
// This class is the base class for Tour.

VertexList::VertexList() = default;

VertexList::VertexList(std::map<unsigned int, std::pair<unsigned int, unsigned int>> neighbors) :
        neighbors(std::move(neighbors)) {}

// Which vertex comes after current, when previous comes before it?
// previous is necessary to determine the direction
unsigned int VertexList::next(unsigned int previous, unsigned int current) const {
    std::pair<unsigned int, unsigned int> currentNeighbors = neighbors.at(current);
    if (previous != currentNeighbors.first and previous != currentNeighbors.second) {
        throw std::runtime_error(
                "Tour::next: previous (" + std::to_string(previous) + ") is not a neighbor of current (" +
                std::to_string(current) + ")");
    } else if (previous == currentNeighbors.first) {
        return currentNeighbors.second;
    } else { // previous == currentNeighbors.second
        return currentNeighbors.first;
    }
}

// Returns some vertex, that is a neighbor of current, if one exists
unsigned int VertexList::next(unsigned int current) const {
    std::pair<unsigned int, unsigned int> currentNeighbors = neighbors.at(current);
    if (currentNeighbors.first != NO_VERTEX) {
        return currentNeighbors.first;
    } else if (currentNeighbors.second != NO_VERTEX) {
        return currentNeighbors.second;
    } else {
        return NO_VERTEX;
    }
}

// Set the vertex that comes after current, when previous comes before it. The neighbors of next are not changed!
// previous is necessary to determine the direction
void VertexList::setNext(unsigned int previous, unsigned int current, unsigned int next) {
    if (previous == next) {
        throw std::runtime_error("Tour::setNext: You cannot set both neighbors to the same vertex.");
    }
    std::pair<unsigned int, unsigned int> currentNeighbors = neighbors.at(current);
    if (previous != currentNeighbors.first and previous != currentNeighbors.second) {
        throw std::runtime_error(
                "Tour::setNext: previous (" + std::to_string(previous) + ") is not a neighbor of current (" +
                std::to_string(current) + ")");
    }
    neighbors.at(current) = std::make_pair(previous, next);
}

// Tries to make vertex1 a neighbor of vertex2 and vertex2 a neighbor of vertex1 and returns wheter this was successful.
// It only fails if one of them already has two neighbors.
bool VertexList::makeNeighbors(unsigned int vertex1, unsigned int vertex2) {
    std::pair<unsigned int, unsigned int> neighbors1 = neighbors.at(vertex1);
    std::pair<unsigned int, unsigned int> neighbors2 = neighbors.at(vertex2);
    if (neighbors1.first == NO_VERTEX and neighbors2.first == NO_VERTEX) {
        neighbors.at(vertex1).first = vertex2;
        neighbors.at(vertex2).first = vertex1;
        return true;
    } else if (neighbors1.first == NO_VERTEX and neighbors2.second == NO_VERTEX) {
        neighbors.at(vertex1).first = vertex2;
        neighbors.at(vertex2).second = vertex1;
        return true;
    } else if (neighbors1.second == NO_VERTEX and neighbors2.first == NO_VERTEX) {
        neighbors.at(vertex1).second = vertex2;
        neighbors.at(vertex2).first = vertex1;
        return true;
    } else if (neighbors1.second == NO_VERTEX and neighbors2.second == NO_VERTEX) {
        neighbors.at(vertex1).second = vertex2;
        neighbors.at(vertex2).second = vertex1;
        return true;
    } else {
        return false;
    }
}

// =================================================== Tour class ======================================================

// This class represents a Tour, but as a subclass of VertexList, neighbors can be manipulated in a variety of ways so
// that it is not actually a class. It is possible to check that with "isHamiltonianTour"

// Initialize the neighbors map with a std::list that represents a tour. For each vertex x the neighbors are the
// adjacent entries in "vertexList" (start and end are also adjancent)
// This expects a list containing the numbers from 0 to tour.size()-1 and clears neighbors
void Tour::setVertices(const std::list<unsigned int> &vertexList) {
    neighbors.clear();
    auto it = vertexList.begin();
    unsigned int previous = *it;
    std::advance(it, 1);
    unsigned int current = *it;
    std::advance(it, 1);
    unsigned int next;
    neighbors.emplace(previous, std::make_pair(vertexList.back(), current));
    for (; it != vertexList.end(); ++it) {
        next = *it;
        neighbors.emplace(current, std::make_pair(previous, next));
        previous = current;
        current = next;
    }
    neighbors.emplace(current, std::make_pair(previous, vertexList.front()));
}

Tour::Tour(const std::list<unsigned int> &vertexList) {
    setVertices(vertexList);
}

// Computes the length of the tour with the cost matrix given by a TSPLIB problem
unsigned int Tour::length(TsplibProblem &tsplibProblem) {
    unsigned int sum = 0;
    TourWalker tourWalker(*this, 0);
    unsigned int currentVertex = tourWalker.getNextVertex();
    unsigned int nextVertex = tourWalker.getNextVertex();
    do {
        sum += tsplibProblem.dist(currentVertex, nextVertex);
        currentVertex = nextVertex;
        nextVertex = tourWalker.getNextVertex();
    } while (currentVertex != 0);
    return sum;
}

// Checks if this Tour really is a hamiltonian tour
bool Tour::isHamiltonianTour() {
    unsigned int dimension = neighbors.size();
    std::vector<bool> visited(dimension, false);
    TourWalker tourWalker(*this, 0);
    unsigned int currentVertex = tourWalker.getNextVertex();
    do {
        if (currentVertex >= visited.size() or visited.at(currentVertex)) {
            return false;
        }
        visited.at(currentVertex) = true;
        currentVertex = tourWalker.getNextVertex();
    } while (currentVertex != 0);
    for (bool v : visited) {
        if (!v) return false;
    }
    return true;
}

// Output the tour to the stream out, typically used with std::cout
std::ostream &operator<<(std::ostream &out, const Tour &tour) {
    std::string output;
    TourWalker tourWalker(tour, 0);
    unsigned int currentVertex = tourWalker.getNextVertex();
    do {
        output += std::to_string(currentVertex + 1) + ", ";
        currentVertex = tourWalker.getNextVertex();
    } while (currentVertex != 0);
    out << output.substr(0, output.length() - 2) << std::endl;
    return out;
}


// =============================================== TourWalker class ====================================================

// This class provides a way to walk through a Tour. It stores the tour, its current and next vertex (to store the
// direction) and can be advanced by "getNextVertex". The tour object can still be altered while walking through it!
// TODO: Check the claim above

// Create a TourWalker that starts the walk at vertex first
TourWalker::TourWalker(const Tour &tour, unsigned int first) : TourWalker(tour, first, tour.next(first)) {}

// Create a TourWalker that starts the walk at vertex first and then walks in the direction of vertex second
TourWalker::TourWalker(const Tour &tour, unsigned int first, unsigned int second) : tour(std::move(tour)),
                                                                                    current(first), next(second) {}

// Advance the walk and get the next vertex
unsigned int TourWalker::getNextVertex() {
    unsigned int previous = current;
    current = next;
    next = tour.next(previous, current);
    return previous;
}


// ============================================ TsplibProblem class ====================================================

// This class reads and stores problems from the TSPLIB library, while also checking for syntax errors, logical errors
// or unsupported keywords

// Interpret a single keyword value pair from the specification part of a TSPLIB file
// Expects keyword and value to not have any superfluous whitespaces or line breaks
// Returns an error message if an error occured and an empty string otherwise
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

// Interpret the file "inputFile" as a TSPLIB file and store the information given there
// Returns an error message if an error occured and an empty string otherwise
std::string TsplibProblem::readFile(std::ifstream &inputFile) {
    std::string line;
    std::string lastDataKeyword;
    std::string::size_type delimiterIndex;
    std::vector<int> numbers; // All numbers in the TSPLIB file as they appear in EDGE_WEIGHT_SECTION

    while (inputFile) {
        getline(inputFile, line);
        line = trim(line);
        if (line.empty()) {
            // ignore empty lines
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
                for (unsigned int i = 0; i < dimension; ++i) {
                    for (unsigned int j = 0; j < dimension; ++j) {
                        matrix[i][j] = numbers.at(numbersIndex);
                        numbersIndex++;
                    }
                }
            } else if (edgeWeightFormat == "LOWER_DIAG_ROW") {
                for (unsigned int i = 0; i < dimension; ++i) {
                    for (unsigned int j = 0; j <= i; ++j) {
                        matrix[i][j] = numbers.at(numbersIndex);
                        matrix[j][i] = numbers.at(numbersIndex);
                        numbersIndex++;
                    }
                }
            } else if (edgeWeightFormat == "UPPER_DIAG_ROW") {
                for (unsigned int i = 0; i < dimension; ++i) {
                    for (unsigned int j = i; j < dimension; ++j) {
                        matrix[i][j] = numbers.at(numbersIndex);
                        matrix[j][i] = numbers.at(numbersIndex);
                        numbersIndex++;
                    }
                }
            } else if (edgeWeightFormat == "UPPER_ROW") {
                // The diagonal is never touched, so it is filled with zeros from the initialization
                for (unsigned int i = 0; i < dimension; ++i) {
                    for (unsigned int j = i + 1; j < dimension; ++j) {
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

// Returns the name of the TSPLIB problem
const std::string &TsplibProblem::getName() const {
    return name;
}

// Returns the type of the TSPLIB problem (always TSP)
const std::string &TsplibProblem::getType() const {
    return type;
}

// Returns the number of vertices in the TSPLIB problem
unsigned int TsplibProblem::getDimension() const {
    return dimension;
}

// Compute the distance of i and j
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


// ============================================== TsplibTour class =====================================================

// Interpret a single keyword value pair from the specification part of a TSPLIB tour file
// Expects keyword and value to not have any superfluous whitespaces or line breaks
// Returns an error message if an error occured and an empty string otherwise
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

// Interpret the file "inputFile" as a TSPLIB tour file and store the information given there
// Returns an error message if an error occured and an empty string otherwise
std::string TsplibTour::readFile(std::ifstream &inputFile) {
    std::string line;
    std::string lastDataKeyword;
    std::string::size_type delimiterIndex;
    std::list<unsigned int> tourList;

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
                            tourList.push_back(n - 1);
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
