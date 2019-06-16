//
// Created by Karl Welzel on 14.05.19.
//

#include <algorithm>
#include <vector>
#include <tuple>
#include <utility>
#include <numeric>
#include "LinKernighanHeuristic.h"
#include "PrimsAlgorithm.h"

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
    return result;
}

CandidateEdges CandidateEdges::nearestNeighbors(const TsplibProblem &problem, size_t k) {
    std::vector<vertex_t> allVertices(problem.getDimension());
    std::iota(allVertices.begin(), allVertices.end(), 0);

    CandidateEdges result;
    result.resize(problem.getDimension());
    for (vertex_t v = 0; v < problem.getDimension(); ++v) {
        result[v].resize(k);
        std::partial_sort_copy(allVertices.begin(), allVertices.end(), result[v].begin(), result[v].end(),
                               [&problem, v](vertex_t w1, vertex_t w2) {
                                   return problem.dist(v, w1) < problem.dist(v, w2);
                               });
    }
    return result;
}

CandidateEdges CandidateEdges::alphaNearestNeighbors(const TsplibProblem &problem, size_t k) {
    dimension_t dimension = problem.getDimension();
    std::vector<vertex_t> allVertices(problem.getDimension());
    std::iota(allVertices.begin(), allVertices.end(), 0);

    std::vector<vertex_t> parent;
    std::vector<vertex_t> topologicalOrder;

    // Generate a minimum spanning tree for vertices 0, ..., dimension-2
    auto distFunction = [&problem](vertex_t i, vertex_t j) { return problem.dist(i, j); };
    std::tie(parent, topologicalOrder) = primsAlgorithm(dimension - 1, distFunction);

    // Find the edges to be added to form a 1-tree
    vertex_t special = dimension - 1; // The special node "1" in the 1-tree
    std::vector<vertex_t> specialNeighbors(2);
    auto specialDistanceComparator = [&problem, special](vertex_t w1, vertex_t w2) {
        return problem.dist(special, w1) < problem.dist(special, w2);
    };
    std::partial_sort_copy(topologicalOrder.begin(), topologicalOrder.end(), specialNeighbors.begin(),
                           specialNeighbors.end(), specialDistanceComparator);

    // Initialize the beta array
    std::vector<std::vector<distance_t>> beta(dimension, std::vector<distance_t>(dimension, 0));

    // Set the beta values for edges (special, v)
    for (vertex_t vertex : topologicalOrder) {
        if (vertex == specialNeighbors[0]) {
            // alpha(special, vertex) = 0 so beta(special, vertex) = dist(special, vertex)
            beta[special][vertex] = beta[vertex][special] = problem.dist(special, vertex);
        } else {
            // When (special, vertex) is required to be in the 1-tree then the edge incident with special with highest
            // distance needs to be removed, namely (special, specialNeighbors[1])
            beta[special][vertex] = beta[vertex][special] = problem.dist(special, specialNeighbors[1]);
        }
    }

    // Compute the beta values for all other edges as described in the paper by Keld Helsgaun from 2000
    for (size_t i = 0; i < topologicalOrder.size(); ++i) {
        for (size_t j = i + 1; j < topologicalOrder.size(); ++j) {
            // beta[i][i] has the smallest possible value 0, so it will never be chosen here
            beta[i][j] = beta[j][i] = std::max(beta[i][parent[j]], problem.dist(j, parent[j]));
        }
    }

/*
    for (vertex_t v = 0; v < dimension; ++v) {
        std::cout << v << " : ";
        for (vertex_t w : adjacentVertices[v]) {
            std::cout << w << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "==========" << std::endl;
    for (vertex_t v : topologicalOrder) {
        std::cout << v << ", parent: " << parent[v] << std::endl;
    }
    std::cout << "==========" << std::endl;
    for (vertex_t v = 0; v < dimension; ++v) {
        for (vertex_t w = 0; w < dimension; ++w) {
            std::cout << tsplibProblem.dist(v, w) - beta[v][w] << ", ";
        }
        std::cout << std::endl;
    }
*/

    // Store the k alpha-nearest neighbors for each vertex
    CandidateEdges result;
    result.resize(problem.getDimension());
    for (vertex_t v = 0; v < problem.getDimension(); ++v) {
        result[v].resize(k);
        auto alphaCompare = [&problem, &beta, v](vertex_t w1, vertex_t w2) {
            distance_t w1Distance = problem.dist(v, w1);
            distance_t w2Distance = problem.dist(v, w2);
            return std::make_tuple(w1Distance - beta[v][w1], w1Distance) <
                   std::make_tuple(w2Distance - beta[v][w2], w2Distance);
        };
        std::partial_sort_copy(allVertices.begin(), allVertices.end(), result[v].begin(), result[v].end(),
                               alphaCompare);
    }

    return result;
}

// ============================================= linKernighanHeuristic =================================================

Tour linKernighanHeuristic(const TsplibProblem &tsplibProblem, const Tour &startTour, CandidateEdges &candidateEdges,
                           size_t backtrackingDepth, size_t infeasibilityDepth) {

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
                vertex_t xiPredecessor = currentTour.predecessor(xi);
                vertex_t xiSuccessor = currentTour.successor(xi);
                for (vertex_t x : candidateEdges[xi]) {
                    if (x != currentWalk.at(0)
                        and x != xiPredecessor and x != xiSuccessor // equivalent to !currentTour.containsEdge(xi, x)
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
                        // currentWalk.appendAndClose(neighbor) is not a valid alternating walk if {neighbor, x_0} is an
                        // edge in currentWalk, but this is only possible if neighbor is x_1, so we only need to exclude
                        // this special case
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

