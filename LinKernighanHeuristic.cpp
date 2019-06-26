//
// Created by Karl Welzel on 14.05.19.
//

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include "Tour.h"
#include "TsplibUtils.h"
#include "LinKernighanHeuristic.h"
#include "AlphaDistances.h"

// ============================================= AlternatingWalk class =================================================

AlternatingWalk AlternatingWalk::close() const {
    AlternatingWalk result(*this);
    result.push_back(at(0));
    return result;
}

AlternatingWalk AlternatingWalk::appendAndClose(vertex_t vertex) const {
    AlternatingWalk result(*this);
    result.push_back(vertex);
    result.push_back(at(0));
    return result;
}

bool AlternatingWalk::containsEdge(vertex_t vertex1, vertex_t vertex2) const {
    for (dimension_t i = 0; i < size() - 1; ++i) {
        if ((operator[](i) == vertex1 and operator[](i + 1) == vertex2) or
            (operator[](i) == vertex2 and operator[](i + 1) == vertex1)) {
            return true;
        }
    }
    return false;
}

std::ostream &operator<<(std::ostream &out, const AlternatingWalk &walk) {
    std::string output;
    for (vertex_t vertex : walk) {
        output += std::to_string(vertex) + ", ";
    }
    out << output.substr(0, output.length() - 2);
    return out;
}


// ============================================= CandidateEdges class =================================================

CandidateEdges CandidateEdges::allNeighbors(const TsplibProblem &problem) {
    CandidateEdges result;
    std::vector<vertex_t> allVertices(problem.getDimension());
    std::iota(allVertices.begin(), allVertices.end(), 0);
    result.assign(problem.getDimension(), allVertices);
    for (vertex_t v : allVertices) {
        // Delete v from its neighbors
        result[v].erase(std::remove(result[v].begin(), result[v].end(), v), result[v].end());
    }
    return result;
}

CandidateEdges CandidateEdges::nearestNeighbors(const TsplibProblem &problem, size_t k) {
    std::vector<vertex_t> allVertices(problem.getDimension());
    std::iota(allVertices.begin(), allVertices.end(), 0);

    CandidateEdges result;
    result.resize(problem.getDimension());
    for (vertex_t v = 0; v < problem.getDimension(); ++v) {
        // Sort the k nearest neighbors of v by distance to v and put them in result[v]
        result[v].resize(k);
        std::partial_sort_copy(allVertices.begin(), std::remove(allVertices.begin(), allVertices.end(), v),
                               result[v].begin(), result[v].end(),
                               [&problem, v](vertex_t w1, vertex_t w2) {
                                   return problem.dist(v, w1) < problem.dist(v, w2);
                               });
    }
    return result;
}

CandidateEdges CandidateEdges::alphaNearestNeighbors(const TsplibProblem &problem, size_t k) {
    dimension_t dimension = problem.getDimension();
    auto dist = [&problem](vertex_t i, vertex_t j) { return problem.dist(i, j); };

    // Compute the alpha distances
    std::vector<std::vector<distance_t>> alpha = alphaDistances(dimension, dist);

    std::vector<vertex_t> allVertices(dimension);
    std::iota(allVertices.begin(), allVertices.end(), 0);

    // Store the k alpha-nearest neighbors for each vertex
    CandidateEdges result;
    result.resize(problem.getDimension());
    for (vertex_t v = 0; v < problem.getDimension(); ++v) {
        result[v].resize(k);
        auto alphaCompare = [&dist, &alpha, v](vertex_t w1, vertex_t w2) {
            return std::make_tuple(alpha[v][w1], dist(v, w1)) < std::make_tuple(alpha[v][w2], dist(v, w2));
        };
        std::partial_sort_copy(allVertices.begin(), std::remove(allVertices.begin(), allVertices.end(), v),
                               result[v].begin(), result[v].end(), alphaCompare);
    }

    return result;
}

// ============================================= linKernighanHeuristic =================================================

LinKernighanHeuristic::LinKernighanHeuristic(TsplibProblem &tsplibProblem, CandidateEdges::Type candidateEdgeType)
        : tsplibProblem(tsplibProblem) {
    switch (candidateEdgeType) {
        case CandidateEdges::ALL_NEIGHBORS:
            candidateEdges = CandidateEdges::allNeighbors(tsplibProblem);
            break;
        case CandidateEdges::NEAREST_NEIGHBORS:
            candidateEdges = CandidateEdges::nearestNeighbors(tsplibProblem);
            break;
        case CandidateEdges::ALPHA_NEAREST_NEIGHBORS:
            candidateEdges = CandidateEdges::alphaNearestNeighbors(tsplibProblem);
            break;
    }
}

vertex_t LinKernighanHeuristic::chooseRandomElement(const std::vector<vertex_t> &elements) {
    std::random_device randomNumberGenerator;
    std::uniform_int_distribution<size_t> distribution(0, elements.size() - 1);
    return elements[distribution(randomNumberGenerator)];
}

Tour LinKernighanHeuristic::generateRandomTour() {
    // Initialize a array with all vertices as the remaining vertices (that are to be placed on the tour)
    std::vector<vertex_t> remainingVertices(tsplibProblem.getDimension());
    std::iota(remainingVertices.begin(), remainingVertices.end(), 0);

    // This variable store the order of the vertices on the tour
    std::vector<vertex_t> tourOrder;

    // Start with a random vertex
    vertex_t currentVertex = chooseRandomElement(remainingVertices);
    remainingVertices.erase(remainingVertices.begin() + currentVertex); // The index of currentVertex is currentVertex
    tourOrder.push_back(currentVertex);


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
        tourOrder.push_back(currentVertex);
    }

    return Tour(tourOrder);
}

Tour LinKernighanHeuristic::improveTour(const Tour &startTour) {
    const dimension_t dimension = tsplibProblem.getDimension();

    // TODO: Introduce the "Don't look" bit

    Tour currentTour = startTour;
    std::vector<std::vector<vertex_t>> vertexChoices;
    AlternatingWalk currentWalk;
    AlternatingWalk bestAlternatingWalk;
    signed_distance_t highestGain = 0;

    while (true) {
        // Create set X_0 with all vertices
        vertexChoices.clear();
        vertexChoices.emplace_back(dimension);
        std::iota(vertexChoices[0].begin(), vertexChoices[0].end(), 0);

        currentWalk.clear();
        bestAlternatingWalk.clear();
        highestGain = 0;
        size_t i = 0;
        while (true) {
            if (vertexChoices[i].empty()) {
                if (highestGain > 0) {
                    //distance_t previousLength = tsplibProblem.length(currentTour);
                    //std::cout << "Exchange done: " << bestAlternatingWalk << std::endl;
                    currentTour.exchange(bestAlternatingWalk);
                    //std::cout << " with gain: " << tsplibProblem.exchangeGain(bestAlternatingWalk)
                    //          << " (highestGain = " << highestGain << ")" << std::endl;
                    //std::cout << " new tour: " << currentTour << std::endl;
                    //std::cout << " previous length: " << previousLength << std::endl;
                    //std::cout << " new length: " << tsplibProblem.length(currentTour) << std::endl;
                    break;
                } else { // highestGain == 0
                    if (i == 0) {
                        return currentTour;
                    } else {
                        i = std::min(i - 1, backtrackingDepth);
                        vertexChoices.erase(vertexChoices.begin() + i + 1, vertexChoices.end());
                        currentWalk.erase(currentWalk.begin() + i, currentWalk.end());
                        continue;
                    }
                }
            }

            currentWalk.push_back(vertexChoices[i].back());
            vertexChoices[i].pop_back();

            // TODO: Remove these debugging checks
            // DEBUG:
            if (vertexChoices.size() != i + 1) {
                throw std::runtime_error(
                        "vertexChoices.size() (=" + std::to_string(vertexChoices.size()) + ") is not i+1 (=" +
                        std::to_string(i + 1) + ")");
            }
            if (currentWalk.size() != i + 1) {
                throw std::runtime_error(
                        "currentWalk.size() (=" + std::to_string(currentWalk.size()) + ") is not i+1 (=" +
                        std::to_string(i + 1) + ")");
            }

            if (i % 2 == 1 and i >= 3) {
                AlternatingWalk closedWalk = currentWalk.close(); // closedWalk = (x_0, x_1, ..., x_i, x_0)
                signed_distance_t gain = tsplibProblem.exchangeGain(closedWalk);
                if (gain > highestGain and currentTour.isTourAfterExchange(closedWalk)) {
                    //std::cout << "New highest gain: " << gain << ", value before: " << highestGain << std::endl;
                    bestAlternatingWalk = closedWalk;
                    highestGain = gain;
                }
            }

            vertexChoices.emplace_back(); // Add set X_{i+1}
            vertex_t xi = currentWalk[i];
            if (i % 2 == 1) { // i is odd
                // Determine possible in-edges
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
                // Determine possible out-edges
                // No out-edge should connect back to currentWalk[0], because at this point currentWalk is not a valid
                // alternating walk (even number of elements) and can never be closed in the future
                if (i == 0 and currentBestTour.getDimension() != 0) {
                    // The first edge to be broken may not be on the currently best solution tour
                    vertex_t x0Predecessor = currentBestTour.predecessor(currentWalk[0]);
                    vertex_t x0Successor = currentBestTour.successor(currentWalk[0]);
                    for (vertex_t neighbor : currentTour.getNeighbors(xi)) {
                        if (neighbor != currentWalk[0] and neighbor != x0Predecessor and neighbor != x0Successor) {
                            vertexChoices[i + 1].push_back(neighbor);
                        }
                    }
                } else if (i <= infeasibilityDepth) {
                    for (vertex_t neighbor : currentTour.getNeighbors(xi)) {
                        if (neighbor != currentWalk[0] and !currentWalk.containsEdge(xi, neighbor)) {
                            vertexChoices[i + 1].push_back(neighbor);
                        }
                    }
                } else {
                    for (vertex_t neighbor : currentTour.getNeighbors(xi)) {
                        // currentWalk.appendAndClose(neighbor) is not a valid alternating walk if {neighbor, x_0} is an
                        // edge in currentWalk, but this is only possible if neighbor is x_1, so we only need to exclude
                        // this special case
                        if (neighbor != currentWalk[0]
                            and !currentWalk.containsEdge(xi, neighbor)
                            and neighbor != currentWalk[1]
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

Tour LinKernighanHeuristic::findBestTour(size_t numberOfTrials) {
    if (numberOfTrials < 1) {
        throw std::runtime_error("The number of trials can not be lower than 1.");
    }
    Tour currentTour = improveTour(generateRandomTour());
    currentBestTour = improveTour(generateRandomTour());
    while (--numberOfTrials > 0) {
        currentTour = improveTour(generateRandomTour());
        std::cout << "Trial result: " << tsplibProblem.length(currentTour) << ", best: "
                  << tsplibProblem.length(currentBestTour) << std::endl;
        if (tsplibProblem.length(currentTour) < tsplibProblem.length(currentBestTour)) {
            currentBestTour = currentTour;
        }
    }
    return currentBestTour;
}
