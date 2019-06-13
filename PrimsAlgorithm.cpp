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

std::tuple<std::vector<std::vector<vertex_t>>, std::vector<vertex_t>>
primsAlgorithm(dimension_t dimension, std::function<distance_t(vertex_t, vertex_t)> dist, vertex_t root) {
    distance_t infiniteDistance = std::numeric_limits<distance_t>::max();

    std::vector<std::vector<vertex_t>> adjacentVertices(dimension);
    std::vector<vertex_t> topologicalOrder;

    std::vector<vertex_t> remainingVertices(dimension);
    std::iota(remainingVertices.begin(), remainingVertices.end(), 0);
    std::vector<std::tuple<vertex_t, vertex_t>> cheapestEdge(dimension);
    std::vector<distance_t> cheapestEdgeCost(dimension, infiniteDistance);

    topologicalOrder.push_back(root);
    remainingVertices.erase(std::find(remainingVertices.begin(), remainingVertices.end(), root));
    for (vertex_t otherVertex : remainingVertices) {
        if (dist(root, otherVertex) < cheapestEdgeCost[otherVertex]) {
            cheapestEdge[otherVertex] = std::make_tuple(root, otherVertex);
            cheapestEdgeCost[otherVertex] = dist(root, otherVertex);
        }
    }


    while (!remainingVertices.empty()) {
        auto currentIterator = std::min_element(remainingVertices.begin(), remainingVertices.end(),
                                                [&cheapestEdgeCost](vertex_t v, vertex_t w) {
                                                    return cheapestEdgeCost[v] < cheapestEdgeCost[w];
                                                });
        vertex_t currentVertex = *currentIterator;
        remainingVertices.erase(currentIterator);

        topologicalOrder.push_back(currentVertex);
        vertex_t v, w;
        std::tie(v, w) = cheapestEdge[currentVertex];
        adjacentVertices[v].push_back(w);
        adjacentVertices[w].push_back(v);

        for (vertex_t otherVertex : remainingVertices) {
            if (dist(currentVertex, otherVertex) < cheapestEdgeCost[otherVertex]) {
                cheapestEdge[otherVertex] = std::make_tuple(currentVertex, otherVertex);
                cheapestEdgeCost[otherVertex] = dist(currentVertex, otherVertex);
            }
        }
    }

    return std::make_tuple(adjacentVertices, topologicalOrder);
}
