//
// Created by Karl Welzel on 29.03.2019.
//

#ifndef LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H
#define LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H

#include <utility>
#include <vector>
#include <limits>
#include <list>
#include <iostream>
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


// ================================================ VertexList class ===================================================

// This class maintains a map that stores a maximum of two neighbors of each vertex. There is no inherent direction
// encoded in these neighbors. This gives additional flexibility when compared to std::list while still having similar
// run times. It is possible to create tours that do not contain every vertex or multiple separate paths or even have a
// vertex x with neighbor y while y does not have x as its neighbor, so extra caution is necessary when manipulating
// these objects.
// This class is the base class for Tour.

class VertexList {
protected:
    std::vector<std::vector<vertex_t>> neighbors;

    // Returns the index of neighbor in neighbors[vertex] if present, otherwise -1
    // Because every vertex has no more than 2 neighbors, it always returns -1, 0 or 1
    int neighborIndex(vertex_t vertex, vertex_t neighbor) const;

    // Set newNeighbor as one of the neighbors of vertex.
    // Throws a runtime_error if both of the neighbors of vertex are already set
    void addNeighbor(vertex_t vertex, vertex_t newNeighbor);

    // Deletes neighbor as a neighbor of vertex
    // Throws a runtime_error if neighbor is not a neighbor of vertex
    void removeNeighbor(vertex_t vertex, vertex_t neighbor);

public:
    VertexList();

    // Create dimension vertices, each with no neighbors
    explicit VertexList(dimension_t dimension);

    // Create VertexList from a given neighbors map
    explicit VertexList(std::vector<std::vector<vertex_t>> neighbors);

    // Returns the number of vertices
    dimension_t getDimension() const;

    // Returns the neighbors of vertex in the tour
    std::vector<vertex_t> getNeighbors(vertex_t vertex);

    // Which vertex comes after current, when previous comes before it?
    // previous is necessary to determine the direction
    // Throws a runtime_error if the vertex cannot be determined
    vertex_t next(vertex_t previous, vertex_t current) const;

    // Returns some vertex, that is a neighbor of current
    // Throws a runtime_error if current has no neighbor
    vertex_t next(vertex_t current) const;

    // Checks whether next is a neighbor of current
    bool isNeighbor(vertex_t vertex, vertex_t neighbor) const;

    // Checks whether the tour contains the edge {vertex1, vertex2}
    // This function is an alias for isNeighbor
    bool containsEdge(vertex_t vertex1, vertex_t vertex2) const;

    // Tries to make vertex1 a neighbor of vertex2 and vertex2 a neighbor of vertex1 and returns whether this was
    // successful. It only fails if one of them already has two neighbors.
    bool makeNeighbors(vertex_t vertex1, vertex_t vertex2);

    std::vector<vertex_t> toTourSequence() const;
};


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

Tour ascendingVerticesHeuristic(const TsplibProblem &tsplibProblem);

#endif //LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H
