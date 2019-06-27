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
#include <unordered_set>
#include "AlphaDistances.h"
#include "PrimsAlgorithm.h"


std::vector<std::vector<distance_t>>
betaValues(dimension_t dimension, const std::function<distance_t(vertex_t, vertex_t)> &dist) {
    std::vector<vertex_t> allVertices(dimension);
    std::iota(allVertices.begin(), allVertices.end(), 0);

    std::vector<vertex_t> parent;
    std::vector<vertex_t> topologicalOrder;

    // Generate a minimum spanning tree in the complete graph
    std::tie(parent, topologicalOrder) = primsAlgorithm(dimension, dist);

    // TODO: Use subgradient optimization

    // Determine the special vertex "1" by selecting the leaf with longest second nearest neighbor distance
    std::unordered_set<vertex_t> leafs(allVertices.begin(), allVertices.end());
    for (vertex_t v : allVertices) {
        leafs.erase(parent[v]);
    }
    auto secondNearestNeighbor = [&allVertices, &dist](vertex_t v) {
        std::vector<vertex_t> twoNearestNeighbors(2);
        std::partial_sort_copy(allVertices.begin(), std::remove(allVertices.begin(), allVertices.end(), v),
                               twoNearestNeighbors.begin(), twoNearestNeighbors.end(),
                               [v, &dist](vertex_t w1, vertex_t w2) {
                                   return dist(v, w1) < dist(v, w2);
                               });
        return twoNearestNeighbors[1];
    };
    vertex_t special = *std::max_element(leafs.begin(), leafs.end(),
                                         [&dist, &secondNearestNeighbor](vertex_t v, vertex_t w) {
                                             return dist(v, secondNearestNeighbor(v)) <
                                                    dist(w, secondNearestNeighbor(w));
                                         });
    vertex_t specialNeighbor = secondNearestNeighbor(special);
    topologicalOrder.erase(std::remove(topologicalOrder.begin(), topologicalOrder.end(), special),
                           topologicalOrder.end());

    // The edge (special, specialNeighbor) is added to the spanning tree to make a 1-tree, so that the two neighbors
    // of special are parent[special] and specialNeighbor with
    //     dist(special, parent[special]) < dist(special, specialNeighbor)

    // Initialize the beta array
    std::vector<std::vector<distance_t>> beta(dimension, std::vector<distance_t>(dimension, 0));

    // Set the beta values for edges (special, v)
    for (vertex_t vertex : topologicalOrder) {
        if (vertex == parent[special]) {
            // alpha(special, vertex) = 0 so beta(special, vertex) = dist(special, vertex)
            beta[special][vertex] = beta[vertex][special] = dist(special, vertex);
        } else {
            // When (special, vertex) is required to be in the 1-tree then the edge incident with special with highest
            // distance needs to be removed, namely (special, specialNeighbor)
            beta[special][vertex] = beta[vertex][special] = dist(special, specialNeighbor);
        }
    }

    // Compute the beta values for all other edges as described in the paper by Keld Helsgaun from 2000
    for (auto i = topologicalOrder.begin(); i != topologicalOrder.end(); ++i) {
        for (auto j = i + 1; j != topologicalOrder.end(); ++j) {
            // beta[x][x] has the smallest possible value 0, so it will never be chosen here
            beta[*i][*j] = beta[*j][*i] = std::max(beta[*i][parent[*j]], dist(*j, parent[*j]));
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