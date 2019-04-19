//
// Created by karl on 19.04.19.
//

#include "TsplibProblem.h"
#include <regex>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>
#include <iostream>
#include <utility>
#include "TsplibTour.h"

// ============================================== TsplibTour class =====================================================

// This class reads and stores tours from the TSPLIB library, while also checking for syntax errors, logical errors
// or unsupported keywords

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
    } else if (keyword == "DIMENSION") {
        dimension = static_cast<unsigned int>(stoi(value));
    } // every unknown keyword is ignored
    return "";
}

// Interpret the file "inputFile" as a TSPLIB tour file and store the information given there
// Returns an error message if an error occurred and an empty string otherwise
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