//
// Created by Karl Welzel on 26.06.19.
//

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <tuple>
#include <vector>
#include <unordered_set>
#include "AlphaDistances.h"
#include "PrimsAlgorithm.h"

signed_distance_t OneTree::length(const std::function<signed_distance_t(vertex_t, vertex_t)> &dist) {
    signed_distance_t result = 0;

    // All edges (v, parent[v])
    for (auto v = topologicalOrder.begin() + 1; v != topologicalOrder.end(); ++v) {
        result += dist(*v, parent[*v]);
    }

    // The edge (special, specialNeighbor)
    result += dist(special, specialNeighbor);

    return result;
}

std::vector<signed_distance_t> OneTree::degrees() {
    std::vector<signed_distance_t> result(parent.size(), 0);

    // All edges (v, parent[v])
    for (auto v = topologicalOrder.begin() + 1; v != topologicalOrder.end(); ++v) {
        result[*v]++;
        result[parent[*v]]++;
    }

    // The edge (special, specialNeighbor)
    result[special]++;
    result[specialNeighbor]++;

    return result;
}

OneTree minimumOneTree(dimension_t dimension, const std::function<signed_distance_t(vertex_t, vertex_t)> &dist) {
    std::unordered_set<vertex_t> allVertices{};
    for (vertex_t v = 0; v < dimension; ++v) {
        allVertices.insert(v);
    }

    std::vector<vertex_t> parent;
    std::vector<vertex_t> topologicalOrder;

    // Generate a minimum spanning tree in the complete graph
    std::tie(parent, topologicalOrder) = primsAlgorithm(dimension, dist);

    // Determine the special vertex "1" by selecting the leaf with longest second nearest neighbor distance
    std::unordered_set<vertex_t> leafs(allVertices);
    for (vertex_t v : allVertices) {
        leafs.erase(parent[v]);
    }
    auto secondNearestNeighbor = [&allVertices, &dist, &parent](vertex_t v) {
        // For every leaf v the parent of v is the nearest neighbor of v, so by excluding v and parent[v] the nearest
        // neighbor that is left is automatically the second nearest
        allVertices.erase(v);
        allVertices.erase(parent[v]);
        vertex_t result = *std::min_element(allVertices.begin(), allVertices.end(),
                                            [v, &dist](vertex_t w1, vertex_t w2) { return dist(v, w1) < dist(v, w2); });
        allVertices.insert(v);
        allVertices.insert(parent[v]);
        return result;
    };
    vertex_t special = *std::max_element(leafs.begin(), leafs.end(),
                                         [&dist, &secondNearestNeighbor](vertex_t v, vertex_t w) {
                                             return dist(v, secondNearestNeighbor(v)) <
                                                    dist(w, secondNearestNeighbor(w));
                                         });

    // Determine the second nearest neighbor of special
    vertex_t specialNeighbor = secondNearestNeighbor(special);

    return OneTree{parent, topologicalOrder, special, specialNeighbor};
}

std::vector<std::vector<distance_t>>
alphaDistances(dimension_t dimension, const std::function<signed_distance_t(vertex_t, vertex_t)> &dist) {
    OneTree tree = minimumOneTree(dimension, dist);

    // Initialize the beta array
    std::vector<std::vector<signed_distance_t>> beta(dimension, std::vector<signed_distance_t>(dimension, 0));

    // The value of beta[i][j] is the length of the edge in the 1-tree that needs to be removed when (i, j) is
    // inserted in the 1-tree. This gives the relationship
    //     alpha[i][j] = dist(i, j) - beta[i][j]

    // Set the beta values for edges (special, v)
    for (vertex_t vertex : tree.topologicalOrder) {
        if (vertex == tree.special) {
            continue;
        } else if (vertex == tree.parent[tree.special]) {
            // alpha(special, vertex) = 0 so beta(special, vertex) = dist(special, vertex)
            beta[tree.special][vertex] = beta[vertex][tree.special] = dist(tree.special, vertex);
        } else {
            // When (special, vertex) is required to be in the 1-tree then the edge incident with special with highest
            // distance needs to be removed, namely (special, specialNeighbor)
            beta[tree.special][vertex] = beta[vertex][tree.special] = dist(tree.special, tree.specialNeighbor);
        }
    }

    // Compute the beta values for all other edges as described in the paper by Keld Helsgaun from 2000
    for (auto i = tree.topologicalOrder.begin(); i != tree.topologicalOrder.end(); ++i) {
        if (*i == tree.special) continue;
        for (auto j = i + 1; j != tree.topologicalOrder.end(); ++j) {
            if (*j == tree.special) continue;
            // beta[x][x] has the smallest possible value 0, so it will never be chosen as the maximum
            beta[*i][*j] = beta[*j][*i] = std::max(beta[*i][tree.parent[*j]], dist(*j, tree.parent[*j]));
        }
    }

    // Initialize the alpha array
    std::vector<std::vector<distance_t>> alpha(dimension, std::vector<distance_t>(dimension, 0));

    // Compute the alpha values
    for (std::size_t i = 0; i < dimension; ++i) {
        for (std::size_t j = 0; j < dimension; ++j) {
            alpha[i][j] = dist(i, j) - beta[i][j];
        }
    }

    return alpha;
}

std::vector<std::vector<distance_t>>
optimizedAlphaDistances(dimension_t dimension, const std::function<signed_distance_t(vertex_t, vertex_t)> &dist) {
    // The penalties of each vertex
    std::vector<signed_distance_t> penalties(dimension, 0);

    // The modified distance function
    auto modifiedDist = [&dist, &penalties](vertex_t v, vertex_t w) {
        return dist(v, w) + penalties[v] + penalties[w];
    };

    // The minimum 1-tree
    OneTree tree = minimumOneTree(dimension, modifiedDist);

    // The lower bound on the size of an optimal tour. This is the objective function we want to maximize
    auto objectiveFunction = [&tree, &modifiedDist, &penalties]() {
        signed_distance_t penaltiesSum = 0;
        for (signed_distance_t p : penalties) {
            penaltiesSum += p;
        }
        return tree.length(modifiedDist) - 2 * penaltiesSum;
    };

    std::size_t stepSize = 1;
    std::size_t periodLength = dimension / 2;
    std::size_t iteration = 0; // A counter for the iterations in the current period
    signed_distance_t currentObjective = objectiveFunction(); // The current value of the objective function
    signed_distance_t previousObjective = currentObjective; // The previous value of the objective function
    // Used to double the step size in the first period until the objective function does not increase
    bool doubleStepSize = true;

    // The subgradient vector is the degree of each vertex minus 2
    std::vector<signed_distance_t> currentSubgradient = tree.degrees();
    std::transform(currentSubgradient.begin(), currentSubgradient.end(), currentSubgradient.begin(),
                   [](signed_distance_t d) { return d - 2; });
    std::vector<signed_distance_t> previousSubgradient = currentSubgradient;

    // Stop the subgradient optimization if the step size the length of the period or the gradient vector is zero
    while (stepSize != 0 and periodLength != 0 and !std::all_of(currentSubgradient.begin(), currentSubgradient.end(),
                                                                [](signed_distance_t d) { return d == 0; })) {
        // Start of the period
        while (iteration++ < periodLength and !std::all_of(currentSubgradient.begin(), currentSubgradient.end(),
                                                           [](signed_distance_t d) { return d == 0; })) {
            // Update the penalties
            for (std::size_t i = 0; i < penalties.size(); ++i) {
                penalties[i] += lround(stepSize * (0.7 * currentSubgradient[i] + 0.3 * previousSubgradient[i]));
            }

            // Update the tree and the subgradient vector
            tree = minimumOneTree(dimension, modifiedDist);
            previousSubgradient = currentSubgradient;
            currentSubgradient = tree.degrees();
            std::transform(currentSubgradient.begin(), currentSubgradient.end(), currentSubgradient.begin(),
                           [](signed_distance_t d) { return d - 2; });
            currentObjective = objectiveFunction();

            // In the first period the step size is doubled until the objective function does not increase
            if (doubleStepSize) {
                if (currentObjective > previousObjective) {
                    stepSize *= 2;
                } else {
                    doubleStepSize = false;
                }
            }

            // If the last iteration in the period leads to an increase of the objective function the period length is
            // doubled
            if (iteration + 1 == periodLength and currentObjective > previousObjective) {
                periodLength *= 2;
            }

            previousObjective = currentObjective;
        }

        // End of the period
        stepSize /= 2;
        periodLength /= 2;
        iteration = 0;
    }

    return alphaDistances(dimension, modifiedDist);
}

