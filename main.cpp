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
#include "SimpleHeuristic.h"


int main(int argc, char *argv[]) {
    // Interpret the first commandline argument as input file and open this file
    std::ifstream problemFile;
    if (argc < 2) {
        std::cerr << "No file supplied." << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "    LinKerninghanAlgorithm.exe tsplib_problem.tsp [optimal_tour.opt.tour]" << std::endl;
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

    // The introduction assignment
    Tour tour = simpleHeuristic(problem);

    // Output the best tour found by the algorithm and compare it to the optimal tour if given
    unsigned int length = tour.length(problem);
    std::cout << "This is the shortest tour found:" << std::endl;
    std::cout << tour;
    std::cout << "It is " << length << " units long." << std::endl << std::endl;

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
        std::cout << "It is " << optimalLength << " units long." << std::endl << std::endl;
        std::cout << "The best tour found by the heuristic is " << (length / (double) optimalLength - 1) * 100
                  << "% above the optimum." << std::endl;
    }

    return 0;
}