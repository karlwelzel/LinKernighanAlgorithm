//
// Created by Karl Welzel on 26.06.19.
//

#ifndef LINKERNINGHANALGORITHM_ALPHADISTANCES_H
#define LINKERNINGHANALGORITHM_ALPHADISTANCES_H


#include <functional>
#include <vector>
#include "Tour.h"

// The alpha distance alpha[i][j] between two vertices in a complete graph is defined as the increase in in length of a
// minimum 1-tree when it is required to contain the edge (i, j). This is a measure for the likelihood that the edge
// (i, j) is contained in an optimum tour (smaller = more likely)

// This struct represents a 1-tree
struct OneTree {
    // parent maps each vertex to its parent in the ordinary minimum spanning tree
    std::vector<vertex_t> parent;

    // topologicalOrder contains all vertices except special in an order such that parent[v] comes before v for every
    // vertex v
    std::vector<vertex_t> topologicalOrder;

    // The special "1" vertex
    vertex_t special;

    // (special, specialNeighbor) is the additional edge added to the tree to form a 1-tree
    vertex_t specialNeighbor;

    signed_distance_t length(const std::function<signed_distance_t(vertex_t, vertex_t)> &dist);

    // Computes the degree of each vertex in the minimum 1-tree
    std::vector<signed_distance_t> degrees();
};

// Computes a minimum 1-tree in the complete graph with dimension vertices and edge weights given by dist. This is done
// by computing a minimum spanning tree in the complete graph, choosing the leaf with the longest second nearest
// neighbor distance as the special vertex and adding the edge to special's second nearest neighbor to the tree. special
// is incident to the edges (special, parent[special]) and (special, specialNeighbor)
// It is guaranteed that dist(special, parent[special]) <= dist(special, specialNeighbor)
OneTree minimumOneTree(dimension_t dimension, const std::function<signed_distance_t(vertex_t, vertex_t)> &dist);

// Computes the alpha distances in the complete graph with dimension vertices and edge weights given by dist
std::vector<std::vector<distance_t>>
alphaDistances(dimension_t dimension, const std::function<signed_distance_t(vertex_t, vertex_t)> &dist);

// Computes the alpha distances in the complete graph with dimension vertices and modified edge weights
// The edge weights are modified by determining a penalty for each vector. The penalty increases or decreases the length
// of all edges incident to this vertex. The penalties are chosen in such a way that the minimum 1-tree with these
// modified distances gives a maximum lower bound on the length of an optimum tour.
std::vector<std::vector<distance_t>>
optimizedAlphaDistances(dimension_t dimension, const std::function<signed_distance_t(vertex_t, vertex_t)> &dist);

#endif //LINKERNINGHANALGORITHM_ALPHADISTANCES_H
