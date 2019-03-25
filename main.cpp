#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>

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

    int dist(unsigned int, unsigned int);
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

int TsplibProblem::dist(unsigned int i, unsigned int j) {
    if (edgeWeightType == "EUC_2D") {
        double d = std::hypot(coordinates.at(i).at(0) - coordinates.at(j).at(0),
                              coordinates.at(i).at(1) - coordinates.at(j).at(0));
        return std::lround(d);
    } else if (edgeWeightType == "CEIL_2D") {
        double d = std::hypot(coordinates.at(i).at(0) - coordinates.at(j).at(0),
                              coordinates.at(i).at(1) - coordinates.at(j).at(0));
        return std::ceil(d);
    }
    return 0;
}


int main() {
    std::ifstream inputFile(R"(D:\Windows 10\Downloads\ch130.tsp)");
    TsplibProblem problem(inputFile);
    inputFile.close();
    std::cout << "This is the distance between node 1 and 2: " << problem.dist(1, 2) << ".";
    return 0;
}