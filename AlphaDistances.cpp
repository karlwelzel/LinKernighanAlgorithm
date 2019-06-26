//
// Created by Karl Welzel on 26.06.19.
//

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <numeric>
#include <tuple>
#include <vector>
#include "AlphaDistances.h"
#include "PrimsAlgorithm.h"

std::vector<std::vector<distance_t>>
betaValues(dimension_t dimension, const std::function<distance_t(vertex_t, vertex_t)> &dist) {
    std::vector<vertex_t> allVertices(dimension);
    std::iota(allVertices.begin(), allVertices.end(), 0);

    std::vector<vertex_t> parent;
    std::vector<vertex_t> topologicalOrder;

    // Generate a minimum spanning tree for vertices 0, ..., dimension-2
    std::tie(parent, topologicalOrder) = primsAlgorithm(dimension - 1, dist);

    // TODO: Don't fix the special node
    // TODO: Use subgradient optimization

    // Find the edges to be added to form a 1-tree
    vertex_t special = dimension - 1; // The special node "1" in the 1-tree
    std::vector<vertex_t> specialNeighbors(2);
    auto specialDistanceComparator = [&dist, special](vertex_t w1, vertex_t w2) {
        return dist(special, w1) < dist(special, w2);
    };
    std::partial_sort_copy(topologicalOrder.begin(), topologicalOrder.end(), specialNeighbors.begin(),
                           specialNeighbors.end(), specialDistanceComparator);

    // Initialize the beta array
    std::vector<std::vector<distance_t>> beta(dimension, std::vector<distance_t>(dimension, 0));

    // Set the beta values for edges (special, v)
    for (vertex_t vertex : topologicalOrder) {
        if (vertex == specialNeighbors[0]) {
            // alpha(special, vertex) = 0 so beta(special, vertex) = dist(special, vertex)
            beta[special][vertex] = beta[vertex][special] = dist(special, vertex);
        } else {
            // When (special, vertex) is required to be in the 1-tree then the edge incident with special with highest
            // distance needs to be removed, namely (special, specialNeighbors[1])
            beta[special][vertex] = beta[vertex][special] = dist(special, specialNeighbors[1]);
        }
    }

    // Compute the beta values for all other edges as described in the paper by Keld Helsgaun from 2000
    for (vertex_t i = 0; i < topologicalOrder.size(); ++i) {
        for (vertex_t j = i + 1; j < topologicalOrder.size(); ++j) {
            // beta[i][i] has the smallest possible value 0, so it will never be chosen here
            beta[i][j] = beta[j][i] = std::max(beta[i][parent[j]], dist(j, parent[j]));
        }
    }

    return beta;
}

std::vector<std::vector<distance_t>>
alphaDistances(dimension_t dimension, const std::function<distance_t(vertex_t, vertex_t)> &dist) {
    std::vector<std::vector<distance_t>> beta = betaValues(dimension, dist);
    std::vector<std::vector<distance_t>> alpha(dimension, std::vector<distance_t>(dimension, 0));

    for (size_t i = 0; i < dimension; ++i) {
        for (size_t j = 0; j < dimension; ++j) {
            alpha[i][j] = dist(i, j) - beta[i][j];
        }
    }
    return alpha;
}