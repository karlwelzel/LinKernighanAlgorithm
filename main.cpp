#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include "Tour.h"
#include "TsplibUtils.h"
#include "SimpleHeuristic.h"
#include "SignedPermutation.h"
#include "LinKernighanHeuristic.h"

// TSPLIB test instances:
// berlin52            EUC_2D
// ch130               EUC_2D (with decimal places)
// d198                EUC_2D (coordinates in scientific notation)
// usa13509            EUC_2D (too big)
// pla7397             CEIL_2D
// brg180              UPPER_ROW
// si175               UPPER_DIAG_ROW
// si1032              UPPER_DIAG_ROW
// fri26               LOWER_DIAG_ROW
// swiss42             FULL_MATRIX


int main(int argc, char *argv[]) {
    std::ifstream problemFile;
    if (argc < 2) {
        std::cerr << "No file supplied." << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "    LinKernighanAlgorithm tsplib_problem.tsp [tsplib_problem.opt.tour]" << std::endl;
        return 1;
    } else {
        problemFile.open(argv[1]);
    }

    // Check whether the problem file could be opened
    if (!problemFile.is_open() or !problemFile.good()) {
        std::cerr << "Could not open the TSPLIB file '" << argv[1] << "'" << std::endl;
        return 1;
    }

    // Interpret the problem file and report any syntax errors, logical errors or unsupported keywords
    TsplibProblem problem;
    std::string errorMessage = problem.readFile(problemFile);
    problemFile.close();

    std::cout << "Opened the " << problem.getName() << " TSPLIB file" << std::endl << std::endl;
    if (!errorMessage.empty()) {
        std::cerr << "The TSPLIB file has an invalid format: " << errorMessage
                  << std::endl;
        return 1;
    }

/*
    // DEBUG: Print out all distances in matrix form
    for (vertex_t i = 0; i < problem.getDimension(); ++i) {
        for (vertex_t j = 0; j < problem.getDimension(); ++j) {
            std::cout.width(10);
            std::cout << problem.dist(i, j);
        }
        std::cout << std::endl;
    }
*/

    // Print the length of the tour 1, ..., n
    distance_t ascendingLength = problem.length(ascendingVerticesHeuristic(problem));
    std::cout << "The tour 1, 2, ..., n has length " << ascendingLength << "." << std::endl << std::endl;

    // Use the heuristic from the introduction assignment to get a tour to start with and optimize it with the
    // Lin-Kernighan-heuristic
    const Tour tour = LinKernighanHeuristic(problem).findBestTour();

    // Output the best tour found by the algorithm and compare it to the optimal tour if given
    distance_t length = problem.length(tour);
    std::cout << "The following tour was the shortest tour found. It is " << length << " units long." << std::endl;
    std::cout << TsplibTour(problem.getName() + ".tour", tour).toTsplibTourFile() << std::endl;

    // Check whether a second command line argument is given
    std::ifstream optimalTourFile;
    if (argc >= 3) {
        // Try to open the file
        optimalTourFile.open(argv[2]);
        if (!optimalTourFile.is_open() or !optimalTourFile.good()) {
            std::cerr << "Could not open the TSPLIB tour file '" << argv[2] << "'" << std::endl;
            return 1;
        }

        // Interpret the tour file and report any syntax errors, logical errors or unsupported keywords
        TsplibTour optimalTour;
        std::string tourErrorMessage = optimalTour.readFile(optimalTourFile);
        optimalTourFile.close();

        if (!tourErrorMessage.empty()) {
            std::cerr << "The TSPLIB tour file has an invalid format: " << tourErrorMessage << std::endl;
            return 1;
        }

        // Check whether the TsplibProblem problem and the TsplibTour optimalTour belong to the same problem
        if (optimalTour.getName().substr(0, problem.getName().size()) != problem.getName()) {
            // optimalTour.getName() does not start with problem.getName()
            std::cerr << "The TSPLIB tour file does not belong to the TSPLIB problem file";
            return 1;
        }

        // Compare the length of the optimal tour to the tour returned by the heuristic
        distance_t optimalLength = problem.length(optimalTour);
        std::cout << "The best tour found by the heuristic is " << ((length / (double) optimalLength) - 1) * 100
                  << "% above the optimum of " << optimalLength << "." << std::endl;
    }


    return 0;
}
