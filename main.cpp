//
// Created by Karl Welzel on 25.03.19.
//

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include "LinKernighanHeuristic.h"
#include "Tour.h"
#include "TsplibUtils.h"

// TODO: Replace all size_t by std::size_t because the latter is the C++ way

int main(int argc, char *argv[]) {
    const std::string helpString = R""(
Usage:
    LinKernighanAlgorithm --help
    LinKernighanAlgorithm tsplib_problem.tsp [options]

Options:
    --number-of-trials=integer
        Sets the maximum number of trials (default: 50)
    --candidate-edges=[ALL|NEAREST|ALPHA_NEAREST|OPT_ALPHA_NEAREST]
        Set the choice of candidate edges for each vertex (default: OPT_ALPHA_NEAREST)
            ALL: all incident edges
            NEAREST: shortest k incident edges
            ALPHA_NEAREST: shortest k incident edges by alpha distance
            OPT_ALPHA_NEAREST: shortest k incident edges by alpha distance after subgradient optimization
    --number-of-candidate-edges=integer
        Set the number of candidate edges k for each vertex. Ignored for candidate-edges=ALL. (default: 5)
    --optimum-tour-length=integer
        Set the optimum tour length for the problem. The output will include the percentage error of the solution.
    --acceptable-error=double
        Set the acceptable error compared to the optimum length in percent. The program will stop when a found tour
        is within the acceptable length range. Ignored if --optimum-value is not given. (default: 1)
    --output-to-file
        Output the best tour found to "tsplib_problem.lk.tour"
    --verbose
        Output debug output, useful to estimate the remaining running time. If not given there is no output until the
        result is computed.

Example:
    LinKernighanAlgorithm si1032.tsp --optimum-tour-length=92650 --acceptable-error=0.1 --output-to-file --verbose

Restrictions:
    The dimension of the TSPLIB problem may not be smaller than 3 or the number of candidate edges plus 1.
)"";

    if (argc < 2) {
        std::cerr << "No TSPLIB file was supplied." << std::endl;
        std::cout << helpString;
        return 1;
    } else if (strcmp(argv[1], "--help") == 0) {
        std::cout << helpString;
    }

    size_t numberOfTrials = 50;
    CandidateEdges::Type candidateEdgeType = CandidateEdges::OPTIMIZED_ALPHA_NEAREST_NEIGHBORS;
    size_t numberOfCandidateEdges = 5;
    distance_t optimumTourLength = 0;
    double acceptableError = 1;
    bool outputToFile = false;
    bool verboseOutput = false;

    std::stringstream stringStream;
    std::string option;
    for (int i = 2; i < argc; ++i) {
        stringStream << argv[i];
        if (!std::getline(stringStream, option, '=')) {
            stringStream.clear();
            continue;
        }
        if (option == "--number-of-trials") {
            stringStream >> numberOfTrials;
        } else if (option == "--candidate-edges") {
            std::string type;
            std::getline(stringStream, option);
            if (type == "ALL") {
                candidateEdgeType = CandidateEdges::ALL_NEIGHBORS;
            } else if (type == "NEAREST") {
                candidateEdgeType = CandidateEdges::NEAREST_NEIGHBORS;
            } else if (type == "ALPHA_NEAREST") {
                candidateEdgeType = CandidateEdges::ALPHA_NEAREST_NEIGHBORS;
            } else if (type == "OPT_ALPHA_NEAREST") {
                candidateEdgeType = CandidateEdges::OPTIMIZED_ALPHA_NEAREST_NEIGHBORS;
            } else {
                std::cerr << "The --candidate-edges type '" << type << "' is not valid" << std::endl;
                std::cout << helpString;
                return 1;
            }
        } else if (option == "--number-of-candidate-edges") {
            stringStream >> numberOfCandidateEdges;
        } else if (option == "--optimum-tour-length") {
            stringStream >> optimumTourLength;
        } else if (option == "--acceptable-error") {
            stringStream >> acceptableError;
        } else if (option == "--output-to-file") {
            outputToFile = true;
        } else if (option == "--verbose") {
            verboseOutput = true;
        } else {
            std::cerr << "An unknown option was given" << std::endl;
            std::cout << helpString;
            return 1;
        }
        if (stringStream.fail()) {
            std::cerr << "One of the options has an invalid format" << std::endl;
            std::cout << helpString;
            return 1;
        }
        stringStream.clear();
    }

    std::ifstream problemFile(argv[1]);

    // Check whether the problem file could be opened
    if (!problemFile.is_open() or !problemFile.good()) {
        std::cerr << "Could not open the TSPLIB file '" << argv[1] << "'" << std::endl;
        return 1;
    }

    // Interpret the problem file and report any syntax errors, logical errors or unsupported keywords
    TsplibProblem problem;
    std::string errorMessage = problem.readFile(problemFile);
    problemFile.close();
    if (!errorMessage.empty()) {
        std::cerr << "The TSPLIB file has an invalid format: " << errorMessage
                  << std::endl;
        return 1;
    }

    if (verboseOutput) std::cout << "Opened the " << problem.getName() << " TSPLIB file" << std::endl;

    CandidateEdges candidateEdges = CandidateEdges::create(problem, candidateEdgeType, numberOfCandidateEdges);

    if (verboseOutput) std::cout << "Computed candidate edges" << std::endl;

    LinKernighanHeuristic heuristic(problem, candidateEdges);
    const Tour tour = heuristic.findBestTour(numberOfTrials, optimumTourLength, acceptableError / 100, verboseOutput);

    // Output the best tour found by the algorithm
    std::string tourName = problem.getName() + ".lk.tour";
    if (outputToFile) {
        std::ofstream outputFile(tourName);
        outputFile << TsplibTour(tourName, tour).toTsplibTourFile() << std::endl;
        outputFile.close();

        if (verboseOutput) std::cout << "Successfully written the tour to '" << tourName << "'" << std::endl;
    } else {
        std::cout << std::endl << TsplibTour(tourName, tour).toTsplibTourFile() << std::endl;
    }

    // Compare the tour length to the length of the optimal tour if given
    distance_t length = problem.length(tour);
    if (optimumTourLength == 0) {
        std::cout << "The shortest tour found is " << length << " units long." << std::endl;
    } else {
        std::cout << "The shortest tour found is " << length << " units long, so "
                  << ((length / (double) optimumTourLength) - 1) * 100 << "% above the optimum of " << optimumTourLength
                  << "." << std::endl;
    }

    return 0;
}
