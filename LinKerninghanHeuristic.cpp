//
// Created by karl on 14.05.19.
//

#include "SimpleHeuristic.h"
#include <algorithm>
#include <vector>
#include <tuple>
#include <utility>
#include <numeric>
#include "LinKerninghanHeuristic.h"


bool containsEdge(const std::vector<vertex_t> &walk, vertex_t vertex1, vertex_t vertex2) {
    for (dimension_t i = 0; i < walk.size() - 1; ++i) {
        if ((walk.at(i) == vertex1 and walk.at(i + 1) == vertex2) or
            (walk.at(i) == vertex2 and walk.at(i + 1) == vertex1)) {
            return true;
        }
    }
    return false;
}

// This is implemented as described in Combinatorial Optimization with p_1 = 5, p_2 = 2 and G = K_n

Tour linKerninghanHeuristic(const TsplibProblem &tsplibProblem, const Tour &startTour) {
    const size_t backtrackingDepth = 5;
    const size_t infeasibilityDepth = 2;

    const dimension_t dimension = tsplibProblem.getDimension();
    Tour currentTour = startTour;
    std::vector<std::vector<vertex_t>> vertexChoices;
    std::vector<vertex_t> chosenVertices;
    std::vector<vertex_t> bestAlternatingWalk;
    distance_t highestGain = 0;

    while (true) {
        // Create set X_0 with all vertices
        vertexChoices.clear();
        vertexChoices.emplace_back(dimension);
        std::iota(vertexChoices.at(0).begin(), vertexChoices.at(0).end(), 0);

        bestAlternatingWalk.clear();
        highestGain = 0;
        size_t i = 0;
        while (true) {
            if (vertexChoices.at(i).empty() and highestGain > 0) {
                currentTour.exchange(bestAlternatingWalk);
                break;
            } else if (vertexChoices.at(i).empty() and highestGain == 0) {
                if (i == 0) {
                    return currentTour;
                } else {
                    i = std::min(i - 1, backtrackingDepth);
                    vertexChoices.erase(vertexChoices.begin() + i, vertexChoices.end());
                    chosenVertices.erase(chosenVertices.begin() + i - 1, chosenVertices.end());
                    continue;
                }
            }

            chosenVertices.push_back(vertexChoices.at(i).back());
            vertexChoices.at(i).pop_back();

            // DEBUG:
            if (vertexChoices.size() == i) {
                throw std::runtime_error(
                        "vertexChoices.size() (=" + std::to_string(vertexChoices.size()) + ") is not i-1 (=" +
                        std::to_string(i - 1) + ")");
            }
            if (chosenVertices.size() == i) {
                throw std::runtime_error(
                        "chosenVertices.size() (=" + std::to_string(chosenVertices.size()) + ") is not i-1 (=" +
                        std::to_string(i - 1) + ")");
            }

            if (i % 2 == 1 and i >= 3) {
                chosenVertices.push_back(chosenVertices.at(0)); // chosenVertices = (x_0, x_1, ..., x_i, x_0)
                distance_t gain;
                if ((gain = tsplibProblem.exchangeGain(chosenVertices)) > highestGain and
                    currentTour.isTourAfterExchange(chosenVertices)) {
                    bestAlternatingWalk = chosenVertices; // TODO: This should be a copy operation, is it?
                    highestGain = gain;
                }
                chosenVertices.pop_back();
            }

            vertexChoices.emplace_back(); // Add set X_{i+1}
            vertex_t xi = chosenVertices.at(i);
            if (i % 2 == 1) { // i is odd
                for (vertex_t x = 0; x < dimension; ++x) {
                    if (x != xi and x != chosenVertices.at(0) and !currentTour.containsEdge(xi, x) and
                        !containsEdge(chosenVertices, xi, x) and
                        tsplibProblem.exchangeGain(chosenVertices) - tsplibProblem.dist(xi, x) > highestGain) {
                        vertexChoices.at(i + 1).push_back(x);
                    }

                }
            } else { // i is even
                if (i <= infeasibilityDepth) {
                    for (vertex_t neighbor : currentTour.getNeighbors(xi)) {
                        if (!containsEdge(chosenVertices, xi, neighbor)) {
                            vertexChoices.at(i + 1).push_back(neighbor);
                        }
                    }
                } else {
                    for (vertex_t neighbor : currentTour.getNeighbors(xi)) {
                        if (!containsEdge(chosenVertices, xi, neighbor)) {
                            chosenVertices.push_back(neighbor);
                            chosenVertices.push_back(chosenVertices.at(0));
                            if (currentTour.isTourAfterExchange(chosenVertices)) {
                                vertexChoices.at(i + 1).push_back(neighbor);
                            }
                            chosenVertices.erase(chosenVertices.end() - 2, chosenVertices.end());
                        }
                    }
                }
            }

            ++i;
        }

    }
}
