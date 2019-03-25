#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

// TODO: Implement remaining functions of TsplibProblem, consider moving TsplibProblem to its own file

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
    int dimension;
    std::string edgeWeightType;
    std::string edgeWeightFormat;

    std::vector<int> numbers;              // All numbers in the TSPLIB file as they appear in EDGE_WEIGHT_SECTION
    std::vector<std::vector<int>> matrix;  // The matrix described in the TSPLIB file

    bool interpretKeyword(std::string, std::string);

    bool fillMatrix();

    bool readProblem(std::ifstream &);

public:
    TsplibProblem(std::ifstream &);

    int cost(int, int);
};

bool TsplibProblem::readProblem(std::ifstream &file) {
    std::string line;
    bool matrixPart = false;
    unsigned int delimiterIndex;

    while (file) {
        getline(file, line);
        line = trim(line);
        if ((delimiterIndex = line.find(delimiter)) != std::string::npos) {
            string keyword = trim(line.substr(0, delimiterIndex));
            string value = trim(line.substr(delimiterIndex + 1, std::string::npos));
            if (!interpretKeyword(keyword, value)) {
                return false;
            }
        } else if (line == "EOF") {
            break;
        } else if (line == "EDGE_WEIGHT_SECTION") {
            matrixPart = true;
        } else if (matrixPart) {
            std::stringstream stream(line);
            int n;
            while (stream >> n) {
                numbers.push_back(n);
            }
        } else {
            return false;
        }
    }

    fillMatrix();

    return true;

    int main() {
        std::cout << "Hello, World!" << std::endl;
        return 0;
    }