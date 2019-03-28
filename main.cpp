#include <utility>
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <algorithm>


// ============================================ TsplibProblem class ====================================================
// consider moving TsplibProblem to its own file

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
    std::vector<int> numbers;              // All numbers in the TSPLIB file as they appear in EDGE_WEIGHT_SECTION
    std::vector<std::vector<int>> matrix;  // The matrix described in the TSPLIB file

    // if EDGE_WEIGHT_TYPE is *_2D
    std::vector<std::vector<double>> coordinates; // Every entry is a 2d coordinate

    //bool fillMatrix();
    bool interpretKeyword(const std::string &, const std::string &);

    bool readProblem(std::ifstream &);

public:
    explicit TsplibProblem(std::ifstream &);

    const std::string &getName() const;

    const std::string &getType() const;

    unsigned int getDimension() const;


    int dist(unsigned int, unsigned int) const;
};

bool TsplibProblem::interpretKeyword(const std::string &keyword, const std::string &value) {
    if (keyword == "NAME") {
        name = value;
    } else if (keyword == "TYPE") {
        type = value;
        if (value != "TSP") {
            return false;
        }

    } else if (keyword == "COMMENT") {
    } else if (keyword == "DIMENSION") {
        dimension = static_cast<unsigned int>(stoi(value));
        if (coordinates.size() != dimension) {
            coordinates.clear();
            while (coordinates.size() < dimension) {
                coordinates.emplace_back(); // Adds an empty vector at the of coordinates
            }
        }
    } else if (keyword == "EDGE_WEIGHT_TYPE") {
        edgeWeightType = value;
        if (edgeWeightType != "EUC_2D" and edgeWeightType != "CEIL_2D") {
            return false;
        }
    } else if (keyword == "EDGE_WEIGHT_FORMAT") {
        edgeWeightFormat = value;
    } else if (keyword == "NODE_COORD_TYPE") {
        nodeCoordType = value;
        if (nodeCoordType != "TWOD_COORDS") {
            return false;
        }
    }
    return true;
}

bool TsplibProblem::readProblem(std::ifstream &inputFile) {
    std::string line;
    bool isMatrixInput = false;
    bool isCoordsInput = false;
    unsigned int delimiterIndex;

    while (inputFile) {
        getline(inputFile, line);
        line = trim(line);
        if ((delimiterIndex = line.find(delimiter)) != std::string::npos) {
            std::string keyword = trim(line.substr(0, delimiterIndex));
            std::string value = trim(line.substr(delimiterIndex + 1, std::string::npos));
            if (!interpretKeyword(keyword, value)) {
                return false;
            }
        } else if (line == "EOF") {
            break;
        } else if (line == "NODE_COORD_SECTION") {
            isCoordsInput = true;
            /* } else if (line == "EDGE_WEIGHT_SECTION") {
                isMatrixInput = true; */
        } else if (isCoordsInput) {
            if (nodeCoordType == "TWOD_COORDS") {
                std::stringstream stream(line);
                unsigned int n;
                double x, y;
                stream >> n >> x >> y;
                coordinates.at(n - 1) = std::vector<double>{x, y};
            } else {
                return false;
            }
            /*} else if (isMatrixInput) {
                std::stringstream stream(line);
                int n;
                while (stream >> n) {
                    numbers.push_back(n);
                } */
        } else {
            return false;
        }
    }

    // Error checking
    if (dimension == 0 or edgeWeightType.empty()) {
        return false;
    }

    if (edgeWeightType == "EUC_2D" or edgeWeightType == "MAX_2D" or edgeWeightType == "MAN_2D"
        or edgeWeightType == "CEIL_2D") {
        if (!isCoordsInput or isMatrixInput) {
            return false;
        }
        for (const std::vector<double> &coord : coordinates) {
            if (coord.size() != 2) {
                return false;
            }
        }
    }
    /* if (edgeWeightType == "EXPLICIT") {
        if (!isMatrixInput or isCoordsInput) {
            return false;
        }
        fillMatrix();
    } */

    return true;
}

TsplibProblem::TsplibProblem(std::ifstream &inputFile) {
    if (!readProblem(inputFile)) {
        std::cerr << "The TSPLIB file has an invalid format or contains unsupported keywords.";
    }
}

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
    }
    return 0;
}

// =================================================== Tour class ======================================================

class Tour {
private:
    std::list<unsigned int> vertices; // The vertices of this tour in the correct order

public:
    explicit Tour(const std::list<unsigned int> &tour) : vertices(tour) {};

    const std::list<unsigned int> &getVertices() const;

    const std::list<unsigned int>::iterator iteratorOf(unsigned int vertex);

    const std::list<unsigned int>::iterator after(std::list<unsigned int>::iterator it);

    const std::list<unsigned int>::iterator before(std::list<unsigned int>::iterator it);

    const unsigned int length(TsplibProblem &tsplibProblem);
};


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
        output += std::to_string(vertex) + ", ";
    }
    out << output.substr(0, output.length() - 2) << std::endl;
    return out;
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
        if (root != -1 and tourParts.getTourPartAt(root).size() == tsplibProblem.getDimension()) {
            return Tour(tourParts.getTourPartAt(tourParts.find(edge.first)));
        }
    }
}

int main(int argc, char *argv[]) {
    std::ifstream inputFile;
    if (argc < 2) {
        std::cout << "No file supplied" << std::endl;
        // return 1;
        inputFile.open(R"(D:\Windows 10\Downloads\berlin52.tsp)");
    } else {
        inputFile.open(argv[1]);
    }

    TsplibProblem problem(inputFile);
    inputFile.close();

    Tour tour = simpleHeuristic(problem);
    std::cout << "This is the shortest route found:" << std::endl;
    std::cout << tour;
    std::cout << "It is " << tour.length(problem) << " units long." << std::endl;

    return 0;
}