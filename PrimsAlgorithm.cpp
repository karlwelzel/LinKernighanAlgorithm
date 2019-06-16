//
// Created by Karl Welzel on 13.06.19.
//

#include <tuple>
#include <vector>
#include <functional>
#include <limits>
#include <numeric>
#include <algorithm>
#include "Tour.h"
#include "PrimsAlgorithm.h"

std::tuple<std::vector<vertex_t>, std::vector<vertex_t>>
primsAlgorithm(dimension_t dimension, const std::function<distance_t(vertex_t, vertex_t)> &dist) {
    distance_t infiniteDistance = std::numeric_limits<distance_t>::max();
    vertex_t noNeighbor = std::numeric_limits<vertex_t>::max();

    // Initialize the parent map and the topological order
    std::vector<vertex_t> parent(dimension);
    std::vector<vertex_t> topologicalOrder;

    // Fill remainingVertices with all vertices
    std::vector<vertex_t> remainingVertices(dimension);
    std::iota(remainingVertices.begin(), remainingVertices.end(), 0);
    std::vector<vertex_t> cheapestNeighborInTree(dimension, noNeighbor);
    std::vector<distance_t> cheapestEdgeCost(dimension, infiniteDistance);

    while (!remainingVertices.empty()) {
        // Get the vertex which can be inserted in the tree with minimal cost and delete it from remainingVertices
        auto currentIterator = std::min_element(remainingVertices.begin(), remainingVertices.end(),
                                                [&cheapestEdgeCost](vertex_t v, vertex_t w) {
                                                    return cheapestEdgeCost[v] < cheapestEdgeCost[w];
                                                });
        vertex_t currentVertex = *currentIterator;
        remainingVertices.erase(currentIterator);

        // Add currentVertex to the tree
        topologicalOrder.push_back(currentVertex);
        vertex_t parentVertex = cheapestNeighborInTree[currentVertex];
        if (parentVertex != noNeighbor) {
            parent[currentVertex] = parentVertex;
        }

        // Update cheapestNeighborInTree
        for (vertex_t otherVertex : remainingVertices) {
            if (dist(currentVertex, otherVertex) < cheapestEdgeCost[otherVertex]) {
                cheapestNeighborInTree[otherVertex] = currentVertex;
                cheapestEdgeCost[otherVertex] = dist(currentVertex, otherVertex);
            }
        }
    }

    return std::make_tuple(parent, topologicalOrder);
}
