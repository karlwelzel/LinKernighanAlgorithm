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


// Prim's algorithm for complete graphs

std::tuple<std::vector<std::vector<vertex_t>>, std::vector<vertex_t>, std::vector<vertex_t>>
primsAlgorithm(dimension_t dimension, const std::function<distance_t(vertex_t, vertex_t)> &dist) {
    distance_t infiniteDistance = std::numeric_limits<distance_t>::max();
    vertex_t noNeighbor = std::numeric_limits<vertex_t>::max();

    std::vector<std::vector<vertex_t>> adjacentVertices(dimension);
    std::vector<vertex_t> topologicalOrder;
    std::vector<vertex_t> parent(dimension);

    std::vector<vertex_t> remainingVertices(dimension);
    std::iota(remainingVertices.begin(), remainingVertices.end(), 0);
    std::vector<vertex_t> cheapestNeighbor(dimension, noNeighbor);
    std::vector<distance_t> cheapestEdgeCost(dimension, infiniteDistance);

    while (!remainingVertices.empty()) {
        auto currentIterator = std::min_element(remainingVertices.begin(), remainingVertices.end(),
                                                [&cheapestEdgeCost](vertex_t v, vertex_t w) {
                                                    return cheapestEdgeCost[v] < cheapestEdgeCost[w];
                                                });
        vertex_t currentVertex = *currentIterator;
        remainingVertices.erase(currentIterator);

        topologicalOrder.push_back(currentVertex);
        vertex_t parentVertex = cheapestNeighbor[currentVertex];
        if (parentVertex != noNeighbor) {
            adjacentVertices[parentVertex].push_back(currentVertex);
            adjacentVertices[currentVertex].push_back(parentVertex);
            parent[currentVertex] = parentVertex;
        }

        for (vertex_t otherVertex : remainingVertices) {
            if (dist(currentVertex, otherVertex) < cheapestEdgeCost[otherVertex]) {
                cheapestNeighbor[otherVertex] = currentVertex;
                cheapestEdgeCost[otherVertex] = dist(currentVertex, otherVertex);
            }
        }
    }

    return std::make_tuple(adjacentVertices, topologicalOrder, parent);
}
