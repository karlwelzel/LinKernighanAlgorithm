//
// Created by Karl Welzel on 29/03/2019.
//

#ifndef LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H
#define LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H

#include <utility>
#include <vector>
#include "Tour.h"
#include "TsplibUtils.h"

/*
 * The assignment is to implement a simple heuristic that sorts all of the edges by length and then tries one-by-one
 * to add them to a subgraph that is always a subgraph of a Hamiltonian tour. This implementation uses an enhanced
 * version of a UnionFind structure called "TourParts" to manage the connected components of this subgraph.
 * The TourParts class is a subclass of VertexList and as such maintains a map that maps each vertex to its neighbors
 * without specifying any direction of these neighbors. Each neighbor represents an outgoing edge of the vertex in the
 * subgraph and every vertex has a maximum of two neighbors to ensure that it is a subgraph of a hamiltonian tour (this
 * is not sufficient, but necessary). Using this map the TourParts class stores the current subgraph. Additionally it
 * uses a UnionFind structures (UnionBySize) to store the connected components in this subgraph. In this case every
 * connected component is path, so TourParts also stores both ends of every path in "pathEnds".
 * At the beginning each vertex is a single connected component (or a path of length one) and then "join(x, y)" tries to
 * connect the paths of x and y by adding the edge (x, y) (this works if and only if both x and y are ends of different
 * paths). This is done for every edge sorted by length until there is only one connected component (= a hamiltonian
 * path). At this point "closeTour" adds the unique edge that connects both ends of the path and returns the resulting
 * tour.
 */


// ============================================== TourParts class ======================================================

class TourParts : public VertexList {
private:
    std::vector<vertex_t> parent;
    std::vector<dimension_t> size;
    std::vector<std::pair<vertex_t, vertex_t>> pathEnds;

public:
    // Initialize this TourParts structure with dimension vertices
    explicit TourParts(dimension_t dimension);

    // Find the parent for vertex x
    vertex_t find(vertex_t x);

    // Try to join the connected components of x and y by adding the edge (x, y)
    // The function does not return whether this was successful
    void join(vertex_t x, vertex_t y);

    // Check whether there is only one connected component. In this case this component is a hamiltonian path and can
    // only be closed (to a tour) by a single edge that connects both ends of the path
    bool isTourClosable();

    // Try to close the tour by adding the unique edge that connects both ends of the only connected component. If there
    // are more than one component (i.e. isTourClosable() == false), this function throws an error, so you should always
    // check isTourClosable() before calling this function
    // After this function has been called the subgraph represented by neighbors is also a tour, so this object should
    // not be used after this function was called
    Tour closeTour();
};

Tour simpleHeuristic(const TsplibProblem &tsplibProblem);

#endif //LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H
