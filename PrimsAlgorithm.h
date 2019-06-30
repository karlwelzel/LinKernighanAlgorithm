//
// Created by Karl Welzel on 13.06.19.
//

#ifndef LINKERNIGHANALGORITHM_PRIMSALGORITHM_H
#define LINKERNIGHANALGORITHM_PRIMSALGORITHM_H

#include <functional>
#include <tuple>
#include <vector>
#include "Tour.h"

// Prim's algorithm for complete graphs

// Given the number of vertices n in a graph and a cost/distance function on each edge of the complete graph K_n the
// function uses Prim's algorithm to calculate a minimum spanning tree in K_n. It returns a map of each vertex to its
// parent and a topological order in which each parent is placed before every child. The running time is O(n^2)
std::tuple<std::vector<vertex_t>, std::vector<vertex_t>>
primsAlgorithm(dimension_t dimension, const std::function<signed_distance_t(vertex_t, vertex_t)> &dist);

#endif //LINKERNIGHANALGORITHM_PRIMSALGORITHM_H
