//
// Created by Karl Welzel on 29/03/2019.
//

// This library contains the code for the introduction assignment

#include <utility>
#include <vector>
#include <algorithm>
#include "SimpleHeuristic.h"

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

// Initialize this TourParts structure with dimension vertices
TourParts::TourParts(const dimension_t dimension) : VertexList(dimension) {
    for (dimension_t i = 0; i < dimension; ++i) {
        parent.push_back(i); // parent[i] = i
        size.push_back(1);   // size[i] = 1
        pathEnds.emplace_back(i, i);
    }
}

// Find the parent for vertex x
vertex_t TourParts::find(const vertex_t x) {
    if (parent.at(x) != x) { // path compression
        parent.at(x) = find(parent.at(x));
    }
    return parent.at(x);
}

// Try to join the connected components of x and y by adding the edge (x, y)
// The function does not return whether this was successful
void TourParts::join(const vertex_t x, const vertex_t y) {
    const vertex_t xr = find(x); // root of x
    const vertex_t yr = find(y); // root of y
    if (xr == yr) {
        return; // A path cannot be joined with itself, this would lead to non-Hamiltonian tours
    }
    std::pair<vertex_t, vertex_t> newPathEnds;
    if (x == pathEnds.at(xr).first and y == pathEnds.at(yr).first) {
        newPathEnds = std::make_pair(pathEnds.at(xr).second, pathEnds.at(yr).second);
    } else if (x == pathEnds.at(xr).first and y == pathEnds.at(yr).second) {
        newPathEnds = std::make_pair(pathEnds.at(xr).second, pathEnds.at(yr).first);
    } else if (x == pathEnds.at(xr).second and y == pathEnds.at(yr).first) {
        newPathEnds = std::make_pair(pathEnds.at(xr).first, pathEnds.at(yr).second);
    } else if (x == pathEnds.at(xr).second and y == pathEnds.at(yr).second) {
        newPathEnds = std::make_pair(pathEnds.at(xr).first, pathEnds.at(yr).first);
    } else {
        return; // x and y cannot be joined, because one of them is not an end of its path
    }

    // Now make x and y neighbors of each other
    if (!makeNeighbors(x, y)) {
        throw std::runtime_error("You cannot join neighbors when one of them already has two neighbors.");
    }

    if (size.at(xr) > size.at(yr)) {
        parent.at(yr) = xr;
        size.at(xr) += size.at(yr);
        pathEnds.at(xr) = newPathEnds;
    } else {
        parent.at(xr) = yr;
        size.at(yr) += size.at(xr);
        pathEnds.at(yr) = newPathEnds;
    }
}

// Check whether there is only one connected component. In this case this component is a hamiltonian path and can only
// be closed (to a tour) by a single edge that connects both ends of the path
bool TourParts::isTourClosable() {
    return size.at(find(0)) == getDimension();
}

// Try to close the tour by adding the unique edge that connects both ends of the only connected component. If there are
// more than one component (i.e. isTourClosable() == false), this function throws an error, so you should always check
// isTourClosable() before calling this function
// After this function has been called the subgraph represented by neighbors is also a tour, so this object should not
// be used after this function was called
Tour TourParts::closeTour() {
    // if there is only one root, then its also the root of 0
    const vertex_t root = find(0);
    if (isTourClosable()) {
        // join the last ends to form a tour
        if (!makeNeighbors(pathEnds.at(root).first, pathEnds.at(root).second)) {
            throw std::runtime_error("You cannot join neighbors when one of them already has two neighbors.");
        }
        return Tour(neighbors);
    } else {
        throw std::runtime_error("The tour cannot be closed while there are multiple connected components");
    }
}


// Compare edges by edge cost/distance given by a TsplibProblem
class EdgeCostComparator {
private:
    TsplibProblem tsplibProblem;

public:
    explicit EdgeCostComparator(const TsplibProblem &tsplibProblem) : tsplibProblem(tsplibProblem) {}

    bool operator()(const std::pair<vertex_t, vertex_t> edge1,
                    const std::pair<vertex_t, vertex_t> edge2) const {
        return (tsplibProblem.dist(edge1.first, edge1.second) < tsplibProblem.dist(edge2.first, edge2.second));
    }
};


Tour simpleHeuristic(const TsplibProblem &tsplibProblem) {
    std::vector<std::pair<vertex_t, vertex_t>> edges;
    for (vertex_t i = 0; i < tsplibProblem.getDimension(); ++i) {
        for (vertex_t j = i + 1; j < tsplibProblem.getDimension(); ++j) {
            edges.emplace_back(i, j);
        }
    }

    const EdgeCostComparator edgeComparator(tsplibProblem);
    // Sort edges by edge cost/distance in ascending order
    std::sort(edges.begin(), edges.end(), edgeComparator);

    TourParts tourParts(tsplibProblem.getDimension());
    for (const std::pair<vertex_t, vertex_t> edge : edges) {
        tourParts.join(edge.first, edge.second);
        if (tourParts.isTourClosable()) {
            return tourParts.closeTour();
        }
    }

    // This line will never be reached, because the graph is a complete graph
    throw std::runtime_error("SimpleHeuristic: Although every edge was checked, no tour has been found.");
}
