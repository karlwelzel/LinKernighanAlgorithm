#include <utility>
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <regex>


// ============================================ TsplibProblem class ====================================================
// consider moving TsplibProblem to its own file

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
    std::vector<std::vector<int>> matrix;  // The matrix described in the TSPLIB file

    // if EDGE_WEIGHT_TYPE is *_2D
    std::vector<std::vector<double>> coordinates; // Every entry is a 2d coordinate

    std::string interpretKeyword(const std::string &keyword, const std::string &value);

public:
    TsplibProblem();

    std::string readFile(std::ifstream &inputFile);

    const std::string &getName() const;

    const std::string &getType() const;

    unsigned int getDimension() const;


    int dist(unsigned int, unsigned int) const;
};

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
    }
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
                    unsigned int n;
                    double x, y;
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

class Tour {
protected:
    std::list<unsigned int> vertices; // The vertices of this tour in the correct order

public:
    Tour();

    explicit Tour(const std::list<unsigned int> &tour);

    const std::list<unsigned int> &getVertices() const;

    const std::list<unsigned int>::iterator iteratorOf(unsigned int vertex);

    const std::list<unsigned int>::iterator after(std::list<unsigned int>::iterator it);

    const std::list<unsigned int>::iterator before(std::list<unsigned int>::iterator it);

    const unsigned int length(TsplibProblem &tsplibProblem);
};

Tour::Tour() = default;

Tour::Tour(const std::list<unsigned int> &tour) : vertices(tour) {}

const std::list<unsigned int> &Tour::getVertices() const {
    return vertices;
}

const std::list<unsigned int>::iterator Tour::iteratorOf(unsigned int vertex) {
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        if (*it == vertex) {
            return it;
        }
    }
}

const std::list<unsigned int>::iterator Tour::after(std::list<unsigned int>::iterator it) {
    if (it == vertices.end()) {
        return vertices.begin();
    } else {
        return std::next(it);
    }
}

const std::list<unsigned int>::iterator Tour::before(std::list<unsigned int>::iterator it) {
    if (it == vertices.begin()) {
        return vertices.end();
    } else {
        return std::prev(it);
    }
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

class TsplibTour : public Tour {
private:
    const char delimiter = ':';

    std::string name;
    std::string type;
    unsigned int dimension = 0;

    std::string interpretKeyword(const std::string &keyword, const std::string &value);

public:
    const std::string &getName() const;

    const std::string &getType() const;

    std::string readFile(std::ifstream &inputFile);
};

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
    }
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

// ============================================== TourParts class ======================================================
class TourParts { // Enhanced UnionFind structure
private:
    std::vector<unsigned int> parent;
    std::vector<unsigned int> rank;
    std::vector<std::list<unsigned int>> tourPart;

public:
    explicit TourParts(unsigned int dimension) {
        for (unsigned int i = 0; i < dimension; ++i) {
            parent.push_back(i); // parent[i] = i
            rank.push_back(0);   // rank[i] = 0
            // The tour only consists of the vertex i
            tourPart.push_back(std::list<unsigned int>{i});
        }
    }

    const std::list<unsigned int> &getTourPartAt(unsigned int i) const {
        return tourPart.at(i);
    }

    unsigned int find(unsigned int x) {
        if (parent.at(x) != x) {
            parent.at(x) = find(parent.at(x));
        }
        return parent.at(x);
    }

    int join(unsigned int x, unsigned int y) {
        // Join the parts of a tour that x and y are part of by adding the edge (x, y)
        unsigned int xr = find(x); // root of x
        unsigned int yr = find(y); // root of y
        if (xr == yr) {
            return -1; // A tour part cannot be joined with itself, this would lead to non-Hamiltonian tours
        }
        if (x == tourPart.at(xr).front() and y == tourPart.at(yr).front()) {
            tourPart.at(xr).reverse();
        } else if (x == tourPart.at(xr).front() and y == tourPart.at(yr).back()) {
            tourPart.at(xr).reverse();
            tourPart.at(yr).reverse();
        } else if (x == tourPart.at(xr).back() and y == tourPart.at(yr).front()) {
        } else if (x == tourPart.at(xr).back() and y == tourPart.at(yr).back()) {
            tourPart.at(yr).reverse();
        } else {
            return -1; // x and y cannot be joined, because one of them is not an end of its tour part
        }
        // The tour parts can now be joined by adding tourPart.at(yr) at the end of tourPart.at(xr) or by adding
        // tourPart.at(xr) at the beginning of tourPart.at(yr)
        unsigned int newRoot;
        if (rank.at(xr) > rank.at(yr)) {
            newRoot = xr;
            parent.at(yr) = xr;
            tourPart.at(xr).splice(tourPart.at(xr).end(), tourPart.at(yr));
        } else {
            newRoot = yr;
            parent.at(xr) = yr;
            tourPart.at(yr).splice(tourPart.at(yr).begin(), tourPart.at(xr));
        }
        if (rank.at(xr) == rank.at(yr)) {
            rank.at(yr)++;
        }
        return newRoot;
    }
};

// Introduction Assignment
class EdgeCostComparator {
private:
    TsplibProblem tsplibProblem;

public:
    explicit EdgeCostComparator(TsplibProblem &tsplibProblem) : tsplibProblem(tsplibProblem) {}

    bool operator()(const std::pair<unsigned int, unsigned int> edge1,
                    const std::pair<unsigned int, unsigned int> edge2) const {
        return (tsplibProblem.dist(edge1.first, edge1.second) < tsplibProblem.dist(edge2.first, edge2.second));
    }
};


Tour simpleHeuristic(TsplibProblem &tsplibProblem) {
    std::vector<std::pair<unsigned int, unsigned int>> edges;
    for (int i = 0; i < tsplibProblem.getDimension(); ++i) {
        for (int j = i + 1; j < tsplibProblem.getDimension(); ++j) {
            edges.emplace_back(i, j);
        }
    }

    EdgeCostComparator edgeComparator(tsplibProblem);
    std::sort(edges.begin(), edges.end(), edgeComparator);

    TourParts tourParts(tsplibProblem.getDimension());
    for (std::pair<unsigned int, unsigned int> edge : edges) {
        int root = tourParts.join(edge.first, edge.second);
        if (root >= 0 and tourParts.getTourPartAt(root).size() == tsplibProblem.getDimension()) {
            return Tour(tourParts.getTourPartAt(root));
        }
    }
}

int main(int argc, char *argv[]) {
    std::ifstream problemFile;
    if (argc < 2) {
        std::cout << "No file supplied" << std::endl;
        return 1;
    } else {
        problemFile.open(argv[1]);
    }

    if (!problemFile.is_open() or !problemFile.good()) {
        std::cerr << "Could not open the TSPLIB file '" << argv[1] << "'" << std::endl;
        return 1;
    }

    TsplibProblem problem;
    std::string errorMessage = problem.readFile(problemFile);
    problemFile.close();

    std::cout << "Opened the " << problem.getName() << " TSPLIB file" << std::endl;

    if (!errorMessage.empty()) {
        std::cerr << "The TSPLIB file has an invalid format: " << errorMessage
                  << std::endl;
        return 1;
    }

    Tour tour = simpleHeuristic(problem);
    unsigned int length = tour.length(problem);
    std::cout << "This is the shortest tour found:" << std::endl;
    std::cout << tour;
    std::cout << "It is " << length << " units long." << std::endl;

    std::ifstream optimalTourFile;
    if (argc >= 3) {
        optimalTourFile.open(argv[2]);
        if (!optimalTourFile.is_open() or !optimalTourFile.good()) {
            std::cerr << "Could not open the TSPLIB tour file '" << argv[2] << "'" << std::endl;
            return 1;
        }

        TsplibTour optimalTour;
        std::string tourErrorMessage = optimalTour.readFile(optimalTourFile);
        optimalTourFile.close();

        if (!tourErrorMessage.empty()) {
            std::cerr << "The TSPLIB tour file has an invalid format: " << tourErrorMessage << std::endl;
            return 1;
        }

        unsigned int optimalLength = optimalTour.length(problem);
        std::cout << "This is the optimal tour:" << std::endl;
        std::cout << optimalTour;
        std::cout << "It is " << optimalLength << " units long." << std::endl;
        std::cout << "The best tour found by the heuristic is " << (length / (double) optimalLength - 1) * 100
                  << "% above the optimum.";
    }

    return 0;
}