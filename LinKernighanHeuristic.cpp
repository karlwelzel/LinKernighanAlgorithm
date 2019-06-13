//
// Created by Karl Welzel on 14.05.19.
//

#include <algorithm>
#include <vector>
#include <tuple>
#include <utility>
#include <numeric>
#include "LinKernighanHeuristic.h"

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

// This is implemented as described in Combinatorial Optimization with p_1 = 5, p_2 = 2 and G = K_n

Tour linKernighanHeuristic(const TsplibProblem &tsplibProblem, const Tour &startTour) {
    const size_t backtrackingDepth = 5;
    const size_t infeasibilityDepth = 2;

    const dimension_t dimension = tsplibProblem.getDimension();
    Tour currentTour = startTour;
    std::vector<std::vector<vertex_t>> vertexChoices;
    AlternatingWalk currentWalk;
    AlternatingWalk bestAlternatingWalk;
    signed_distance_t highestGain = 0;

    while (true) {
        // Create set X_0 with all vertices
        vertexChoices.clear();
        vertexChoices.emplace_back(dimension);
        std::iota(vertexChoices.at(0).begin(), vertexChoices.at(0).end(), 0);

        currentWalk.clear();
        bestAlternatingWalk.clear();
        highestGain = 0;
        size_t i = 0;
        while (true) {
            if (vertexChoices.at(i).empty()) {
                if (highestGain > 0) {
                    distance_t previousLength = tsplibProblem.length(currentTour);
                    std::cout << "Exchange done: " << bestAlternatingWalk << std::endl;
                    currentTour.exchange(bestAlternatingWalk);
                    std::cout << " with gain: " << tsplibProblem.exchangeGain(bestAlternatingWalk) << " (highestGain = "
                              << highestGain << ")" << std::endl;
                    std::cout << " new tour: " << currentTour << std::endl;
                    std::cout << " previous length: " << previousLength << std::endl;
                    std::cout << " new length: " << tsplibProblem.length(currentTour) << std::endl;
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

            currentWalk.push_back(vertexChoices.at(i).back());
            vertexChoices.at(i).pop_back();

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
                    std::cout << "New highest gain: " << gain << ", value before: " << highestGain << std::endl;
                    bestAlternatingWalk = closedWalk;
                    highestGain = gain;
                }
            }

            vertexChoices.emplace_back(); // Add set X_{i+1}
            vertex_t xi = currentWalk.at(i);
            if (i % 2 == 1) { // i is odd
                // Determine possible in-edges
                signed_distance_t currentGain = tsplibProblem.exchangeGain(currentWalk);
                // TODO: curentTour.containsEdge, repeatedly calls predecessor and successor of xi
                for (vertex_t x = 0; x < dimension; ++x) {
                    if (x != xi and x != currentWalk.at(0)
                        and !currentTour.containsEdge(xi, x)
                        and !currentWalk.containsEdge(xi, x)
                        and currentGain - static_cast<signed_distance_t>(tsplibProblem.dist(xi, x)) > highestGain) {

                        vertexChoices.at(i + 1).push_back(x);
                    }

                }
            } else { // i is even
                // Determine possible out-edges
                // No out-edge should connect back to currentWalk[0], because at this point currentWalk is not a valid
                // alternating walk (even number of elements) and can never be closed in the future
                if (i <= infeasibilityDepth) {
                    for (vertex_t neighbor : currentTour.getNeighbors(xi)) {
                        if (neighbor != currentWalk.at(0) and !currentWalk.containsEdge(xi, neighbor)) {
                            vertexChoices.at(i + 1).push_back(neighbor);
                        }
                    }
                } else {
                    for (vertex_t neighbor : currentTour.getNeighbors(xi)) {
                        // currentWalk.appendAndClose(neighbor) is not a valid alternating walk if {neighbor, x_0}, but
                        // this is only possible if neighbor is x_1, so we only need to exclude this special case
                        if (neighbor != currentWalk.at(0)
                            and !currentWalk.containsEdge(xi, neighbor)
                            and neighbor != currentWalk.at(1)
                            and currentTour.isTourAfterExchange(currentWalk.appendAndClose(neighbor))) {
                            vertexChoices.at(i + 1).push_back(neighbor);
                        }
                    }
                }
            }

            ++i;
        }

    }
}
