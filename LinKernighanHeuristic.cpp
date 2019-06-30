//
// Created by Karl Welzel on 14.05.19.
//

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>
#include "Tour.h"
#include "TsplibUtils.h"
#include "LinKernighanHeuristic.h"
#include "AlphaDistances.h"

// ============================================= CandidateEdges class =================================================

CandidateEdges::CandidateEdges(dimension_t dimension, const std::vector<vertex_t> &fillValue) : neighbors(dimension,
                                                                                                          fillValue) {}

std::vector<vertex_t> &CandidateEdges::operator[](std::size_t index) {
    return neighbors[index];
}

CandidateEdges CandidateEdges::allNeighbors(const TsplibProblem &problem) {
    std::vector<vertex_t> allVertices(problem.getDimension());
    std::iota(allVertices.begin(), allVertices.end(), 0);

    CandidateEdges result(problem.getDimension(), allVertices);
    for (vertex_t v : allVertices) {
        // Delete v from its neighbors
        result[v].erase(std::remove(result[v].begin(), result[v].end(), v), result[v].end());
    }
    return result;
}

CandidateEdges CandidateEdges::rawNearestNeighbors(
        dimension_t dimension, std::size_t k, const std::function<bool(vertex_t, vertex_t, vertex_t)> &distCompare) {
    std::unordered_set<vertex_t> allVertices{};
    for (vertex_t v = 0; v < dimension; ++v) {
        allVertices.insert(v);
    }

    CandidateEdges result(dimension, std::vector<vertex_t>(k));
    for (vertex_t v = 0; v < dimension; ++v) {
        // Sort the k nearest neighbors of v by distance to v and put them in result[v]
        allVertices.erase(v); // Don't include v in the search
        // std::bind(distCompare, v, std::placeholders::_1, std::placeholders::_2) is a function that takes two vertices
        // and compares them by their distance to v
        std::partial_sort_copy(allVertices.begin(), allVertices.end(), result[v].begin(), result[v].end(),
                               std::bind(distCompare, v, std::placeholders::_1, std::placeholders::_2));
        allVertices.insert(v);
    }
    return result;
}

CandidateEdges CandidateEdges::nearestNeighbors(const TsplibProblem &problem, std::size_t k) {
    auto distCompare = [&problem](vertex_t v, vertex_t w1, vertex_t w2) {
        return problem.dist(v, w1) < problem.dist(v, w2);
    };

    return rawNearestNeighbors(problem.getDimension(), k, distCompare);
}

CandidateEdges CandidateEdges::alphaNearestNeighbors(const TsplibProblem &problem, std::size_t k) {
    auto dist = [&problem](vertex_t i, vertex_t j) { return problem.dist(i, j); };

    // Compute the alpha distances
    std::vector<std::vector<distance_t>> alpha = alphaDistances(problem.getDimension(), dist);

    auto distCompare = [&problem, &alpha](vertex_t v, vertex_t w1, vertex_t w2) {
        return std::make_tuple(alpha[v][w1], problem.dist(v, w1)) <
               std::make_tuple(alpha[v][w2], problem.dist(v, w2));
    };
    return rawNearestNeighbors(problem.getDimension(), k, distCompare);
}

CandidateEdges CandidateEdges::optimizedAlphaNearestNeighbors(const TsplibProblem &problem, std::size_t k) {
    auto dist = [&problem](vertex_t i, vertex_t j) { return problem.dist(i, j); };

    // Compute the optimized alpha distances
    std::vector<std::vector<distance_t>> alpha = optimizedAlphaDistances(problem.getDimension(), dist);

    auto distCompare = [&problem, &alpha](vertex_t v, vertex_t w1, vertex_t w2) {
        return std::make_tuple(alpha[v][w1], problem.dist(v, w1)) <
               std::make_tuple(alpha[v][w2], problem.dist(v, w2));
    };
    return rawNearestNeighbors(problem.getDimension(), k, distCompare);
}

CandidateEdges CandidateEdges::create(const TsplibProblem &problem, CandidateEdges::Type candidateEdgeType,
                                      std::size_t k) {
    switch (candidateEdgeType) {
        case Type::ALL_NEIGHBORS:
            return CandidateEdges::allNeighbors(problem);
        case Type::NEAREST_NEIGHBORS:
            return CandidateEdges::nearestNeighbors(problem, k);
        case Type::ALPHA_NEAREST_NEIGHBORS:
            return CandidateEdges::alphaNearestNeighbors(problem, k);
        default:
        case Type::OPTIMIZED_ALPHA_NEAREST_NEIGHBORS:
            return CandidateEdges::optimizedAlphaNearestNeighbors(problem, k);
    }
}


// ========================================== LinKernighanHeuristic class ==============================================

LinKernighanHeuristic::LinKernighanHeuristic(TsplibProblem &tsplibProblem, CandidateEdges candidateEdges)
        : tsplibProblem(tsplibProblem), candidateEdges(std::move(candidateEdges)) {
}

vertex_t LinKernighanHeuristic::chooseRandomElement(const std::vector<vertex_t> &elements) {
    std::random_device randomNumberGenerator;
    std::uniform_int_distribution<std::size_t> distribution(0, elements.size() - 1);
    return elements[distribution(randomNumberGenerator)];
}

Tour LinKernighanHeuristic::generateRandomTour() {
    // Initialize a array with all vertices as the remaining vertices (that are to be placed on the tour)
    std::vector<vertex_t> remainingVertices(tsplibProblem.getDimension());
    std::iota(remainingVertices.begin(), remainingVertices.end(), 0);

    // This variable stores the order of the vertices on the tour
    std::vector<vertex_t> tourSequence;

    // Start with a random vertex
    vertex_t currentVertex = chooseRandomElement(remainingVertices);
    remainingVertices.erase(remainingVertices.begin() + currentVertex); // The index of currentVertex is currentVertex
    tourSequence.push_back(currentVertex);

    // In each step, decide for each vertex otherVertex if it is an element in one or more of these categories
    // (1) otherVertex was not already chosen and {currentVertex, otherVertex} is a candidate edge and on the
    //     current best tour
    // (2) otherVertex was not already chosen and {currentVertex, otherVertex} is a candidate edge
    // (3) otherVertex was not already chosen
    // Choose the next current vertex randomly from the first non-empty category above

    std::vector<vertex_t> candidatesInBestTour; // Category (1)
    std::vector<vertex_t> candidates; // Category (2)
    // Category (3) is remainingVertices
    while (!remainingVertices.empty()) {
        candidatesInBestTour.clear();
        candidates.clear();
        for (vertex_t otherVertex : candidateEdges[currentVertex]) {
            if (std::find(remainingVertices.begin(), remainingVertices.end(), otherVertex) != remainingVertices.end()) {
                // otherVertex is an element of remainingVertices
                if (currentBestTour.getDimension() != 0 and currentBestTour.containsEdge(currentVertex, otherVertex)) {
                    candidatesInBestTour.push_back(otherVertex);
                }
                candidates.push_back(otherVertex);
            }
        }

        if (!candidatesInBestTour.empty()) {
            currentVertex = chooseRandomElement(candidatesInBestTour);
        } else if (!candidates.empty()) {
            currentVertex = chooseRandomElement(candidates);
        } else {
            currentVertex = chooseRandomElement(remainingVertices);
        }
        remainingVertices.erase(std::remove(remainingVertices.begin(), remainingVertices.end(), currentVertex),
                                remainingVertices.end());
        tourSequence.push_back(currentVertex);
    }

    return Tour(tourSequence);
}

Tour LinKernighanHeuristic::improveTour(const Tour &startTour) {
    const dimension_t dimension = tsplibProblem.getDimension();

    Tour currentTour = startTour;
    // vertexChoices[i] stores all possible choices for vertex x_i. This is used for backtracking
    std::vector<std::vector<vertex_t>> vertexChoices;
    AlternatingWalk currentWalk; // The i-th element of currentWalk is also reffered to as x_i
    AlternatingWalk bestAlternatingWalk;
    signed_distance_t highestGain = 0;

    while (true) {
        // Reset everything
        vertexChoices.clear();
        currentWalk.clear();
        bestAlternatingWalk.clear();
        highestGain = 0;
        std::size_t i = 0;

        // Fill vertexChoices[0] with all vertices
        vertexChoices.emplace_back(dimension);
        std::iota(vertexChoices[0].begin(), vertexChoices[0].end(), 0);

        while (true) {
            if (vertexChoices[i].empty()) {
                // The current alternating walk cannot be expanded further
                if (highestGain > 0) {
                    currentTour.exchange(bestAlternatingWalk);
                    break;
                } else { // highestGain == 0
                    if (i == 0) {
                        // No improvement for currentTour was found
                        return currentTour;
                    } else {
                        // Reset the search to level min(i-1, backtrackingDepth)
                        i = std::min(i - 1, backtrackingDepth);
                        vertexChoices.erase(vertexChoices.begin() + i + 1, vertexChoices.end());
                        currentWalk.erase(currentWalk.begin() + i, currentWalk.end());
                        continue;
                    }
                }
            }

            // Choose x_i from vertexChoices[i] and remove it from there
            vertex_t xi = vertexChoices[i].back();
            vertexChoices[i].pop_back();
            currentWalk.push_back(xi);

            if (i % 2 == 1 and i >= 3) {
                // Check if the exchange of the current walk after closing it leads to a tour and update
                // bestAlternatingWalk if necessary
                AlternatingWalk closedWalk = currentWalk.close(); // closedWalk = (x_0, x_1, ..., x_i, x_0)
                signed_distance_t gain = tsplibProblem.exchangeGain(closedWalk);
                if (gain > highestGain and currentTour.isTourAfterExchange(closedWalk)) {
                    bestAlternatingWalk = closedWalk;
                    highestGain = gain;
                }
            }

            // Add vertexChoices[i+1] and fill it afterwards
            vertexChoices.emplace_back();
            if (i % 2 == 1) { // i is odd
                // Determine possible in-edges (xi, x)
                signed_distance_t currentGain = tsplibProblem.exchangeGain(currentWalk);
                vertex_t xiPredecessor = currentTour.predecessor(xi);
                vertex_t xiSuccessor = currentTour.successor(xi);
                for (vertex_t x : candidateEdges[xi]) {
                    if (x != currentWalk[0]
                        and x != xiPredecessor and x != xiSuccessor // equivalent to !currentTour.containsEdge(xi, x)
                        and !currentWalk.containsEdge(xi, x)
                        and currentGain - static_cast<signed_distance_t>(tsplibProblem.dist(xi, x)) > highestGain) {

                        vertexChoices[i + 1].push_back(x);
                    }
                }
            } else { // i is even
                // Determine possible out-edges (xi, neighbor)

                // For i > infeasibilityDepth the out-edges must be chosen in way that the exchange of the current walk
                // after appending the edge and closing the walk (currentWalk.appendAndClose(neighbor)) leads to tour

                // Special caution is needed because:
                // (1) No out-edge should connect back to x_0, because at this point currentWalk is not a valid
                //     alternating walk (even number of elements) and can never be closed in the future
                // (2) currentWalk.appendAndClose(neighbor) is not a valid alternating walk if {neighbor, x_0} is an
                //     edge in currentWalk, but this is only possible if neighbor = x_1, so we only need to exclude this
                //     special case
                if (i == 0 and currentBestTour.getDimension() != 0) {
                    // The first edge to be broken may not be on the currently best solution tour
                    // (1) can not happen because x_0 is not a neighbor of x_0
                    for (vertex_t neighbor : currentTour.getNeighbors(xi)) {
                        if (!currentBestTour.containsEdge(xi, neighbor)) {
                            vertexChoices[i + 1].push_back(neighbor);
                        }
                    }
                } else if (i <= infeasibilityDepth) {
                    for (vertex_t neighbor : currentTour.getNeighbors(xi)) {
                        if (neighbor != currentWalk[0] // (1)
                            and !currentWalk.containsEdge(xi, neighbor)) {

                            vertexChoices[i + 1].push_back(neighbor);
                        }
                    }
                } else {
                    for (vertex_t neighbor : currentTour.getNeighbors(xi)) {
                        if (neighbor != currentWalk[0] // (1)
                            and !currentWalk.containsEdge(xi, neighbor)
                            and neighbor != currentWalk[1] // (2)
                            and currentTour.isTourAfterExchange(currentWalk.appendAndClose(neighbor))) {

                            vertexChoices[i + 1].push_back(neighbor);
                        }
                    }
                }
            }

            ++i;
        }

    }
}

Tour
LinKernighanHeuristic::findBestTour(std::size_t numberOfTrials, distance_t optimumTourLength, double acceptableError,
                                    bool verboseOutput) {
    if (numberOfTrials < 1) {
        throw std::runtime_error("The number of trials can not be lower than 1.");
    }

    Tour startTour;
    Tour currentTour;
    distance_t currentBestLength = std::numeric_limits<distance_t>::max();
    std::size_t trialCount = 0;

    while (trialCount++ < numberOfTrials) {
        if (verboseOutput) std::cout << "Trial " << trialCount << " | " << std::flush;

        startTour = generateRandomTour();
        if (verboseOutput)
            std::cout << "Length of startTour: " << tsplibProblem.length(startTour) << " | " << std::flush;

        currentTour = improveTour(startTour);
        if (verboseOutput)
            std::cout << "Length of currentTour: " << tsplibProblem.length(currentTour) << " | " << std::flush;

        // Update currentBestTour if necessary
        if (tsplibProblem.length(currentTour) < currentBestLength) {
            currentBestTour = currentTour;
            currentBestLength = tsplibProblem.length(currentBestTour);
        }
        if (verboseOutput) std::cout << "Length of currentBestTour: " << currentBestLength << std::endl;

        // Stop if the increase in length of the current best tour relative to the optimal length is below the
        // threshold set by acceptableError
        if (currentBestLength < (1 + acceptableError) * optimumTourLength) {
            break;
        }
    }

    return currentBestTour;
}
