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

Tour linKerninghanHeuristic(const TsplibProblem &tsplibProblem, const Tour &startTour) {
    const size_t backtrackingDepth = 5;
    const size_t infeasibilityDepth = 2;

    Tour currentTour = startTour;
    std::vector<std::vector<vertex_t>> vertexChoices;
    std::vector<vertex_t> chosenVertices;
    std::vector<vertex_t> bestAlternatingWalk;
    distance_t highestGain = 0;

    while (true) {
        // Create set X_0 with all vertices
        vertexChoices.clear();
        vertexChoices.emplace_back(tsplibProblem.getDimension());
        std::iota(vertexChoices.at(0).begin(), vertexChoices.at(0).end(), 0);

        bestAlternatingWalk.clear();
        highestGain = 0;
        for (size_t i = 0;; ++i) {
            if (vertexChoices.at(i).empty() and highestGain > 0) {
                currentTour.exchange(bestAlternatingWalk);
                break;
            } else if (vertexChoices.at(i).empty() and highestGain == 0) {
                if (i == 0) {
                    return currentTour;
                } else {
                    i = std::min(i - 1, backtrackingDepth);
                    // TODO: Remove unnecessary vertexChoices and chosenVertices
                }
            }

            // TODO: Maybe always remove the last element of chosenVertices to make this if-statement unnecessary
            if (chosenVertices.size() < i) {
                chosenVertices.push_back(vertexChoices.at(i).back());
            } else {
                chosenVertices.at(i) = vertexChoices.at(i).back();
            }
            vertexChoices.at(i).pop_back();

            // DEBUG:
            if (chosenVertices.size() < i) {
                throw std::runtime_error("chosenVertices.size() is too low");
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

            // DEBUG:
            if (vertexChoices.size() != i - 1) {
                throw std::runtime_error("vertexChoices.size() is not correct");
            }
            if (i % 2 == 1) { // i is odd
                vertexChoices.emplace_back(); // Add set X_{i+1}
                // You can use currentTour.isNeighbor()
            } else { // i is even

            }
        }

    }
}