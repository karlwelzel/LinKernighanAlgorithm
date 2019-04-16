//
// Created by Karl Welzel on 29/03/2019.
//

// This library contains the code for the introduction assignment

#include <utility>
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <regex>
#include "TsplibUtils.h"
#include "SimpleHeuristic.h"

// ============================================== TourParts class ======================================================
// Enhanced UnionFind structure

// This class is a VertexList, so it also maintains a map that maps each vertex to its neighbors without specifying any
// direction of these neighbors. At the beginning each vertex has no neighbors and then the individual connected
// components (paths) are connected as needed. To make this as efficient as possible the parent of every connected
// component starts the end vertices of the path that it (and all its children) are in.

TourParts::TourParts(unsigned int dimension) : dimension(dimension) {
    for (unsigned int i = 0; i < dimension; ++i) {
        parent.push_back(i); // parent[i] = i
        size.push_back(1);   // size[i] = 1
        // The tour initially only consists of the vertex i
        neighbors.emplace(i, std::make_pair(VertexList::NO_VERTEX, VertexList::NO_VERTEX));
        pathEnds.emplace_back(i, i);
    }
}

unsigned int TourParts::find(unsigned int x) {
    if (parent.at(x) != x) { // path compression
        parent.at(x) = find(parent.at(x));
    }
    return parent.at(x);
}

void TourParts::join(unsigned int x, unsigned int y) {
    // Join the parts of a tour that x and y are part of by adding the edge (x, y)
    unsigned int xr = find(x); // root of x
    unsigned int yr = find(y); // root of y
    if (xr == yr) {
        return; // A path cannot be joined with itself, this would lead to non-Hamiltonian tours
    }
    std::pair<unsigned int, unsigned int> newPathEnds;
    if (x == pathEnds.at(xr).first and y == pathEnds.at(yr).first) {
        newPathEnds = std::make_pair(pathEnds.at(xr).second, pathEnds.at(yr).second);
    } else if (x == pathEnds.at(xr).first and y == pathEnds.at(yr).second) {
        newPathEnds = std::make_pair(pathEnds.at(xr).second, pathEnds.at(yr).first);
    } else if (x == pathEnds.at(xr).second and y == pathEnds.at(yr).first) {
        newPathEnds = std::make_pair(pathEnds.at(xr).first, pathEnds.at(yr).second);
    } else if (x == pathEnds.at(xr).second and y == pathEnds.at(yr).second) {
        newPathEnds = std::make_pair(pathEnds.at(xr).first, pathEnds.at(yr).first);
    } else {
        return; // x and y cannot be joined, because one of them is not an end of its tour part
    }

    // Now make x and y neighbors of each other
    if (!addNeighbor(x, y) or !addNeighbor(y, x)) {
        throw std::runtime_error("You cannot join neighbors where one of them already has two neighbors.");
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

bool TourParts::isTourClosable() {
    return size.at(find(0)) == dimension;
}

Tour TourParts::closeTour() {
    // if there is only one root, then its also the root of 0
    unsigned int root = find(0);
    if (isTourClosable()) {
        // join the last ends to form a tour
        if (!addNeighbor(pathEnds.at(root).first, pathEnds.at(root).second) or
            !addNeighbor(pathEnds.at(root).second, pathEnds.at(root).first)) {
            throw std::runtime_error("You cannot join neighbors where one of them already has two neighbors.");
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
    explicit EdgeCostComparator(TsplibProblem &tsplibProblem) : tsplibProblem(tsplibProblem) {}

    bool operator()(const std::pair<unsigned int, unsigned int> edge1,
                    const std::pair<unsigned int, unsigned int> edge2) const {
        return (tsplibProblem.dist(edge1.first, edge1.second) < tsplibProblem.dist(edge2.first, edge2.second));
    }
};


/*
 * The assignment was to implement a simple heuristic that sorts all of the edges by length and then tries one-by-one
 * to add them to a subgraph that is always a subgraph of a Hamiltonian tour. This implementation uses an enhanced
 * version of a UnionFind structure called "TourParts" to manage the connected components. In this case every connected
 * component is a path and a part of the final tour, so tourPart stores the current path as a std::list for every head
 * of a component. The "join" function tries to join two vertices and returns whether this was possible (both vertices
 * need to be at the end of a different path). In this way the algorithm tries to join every edge and stops when the
 * new path (after successful joining) contains every vertex. In this case only one edge can be used to close up the
 * tour and the result is returned.
 */
Tour simpleHeuristic(TsplibProblem &tsplibProblem) {
    std::vector<std::pair<unsigned int, unsigned int>> edges;
    for (int i = 0; i < tsplibProblem.getDimension(); ++i) {
        for (int j = i + 1; j < tsplibProblem.getDimension(); ++j) {
            edges.emplace_back(i, j);
        }
    }

    EdgeCostComparator edgeComparator(tsplibProblem);
    // Sort edges by edge cost/distance in ascending order
    std::sort(edges.begin(), edges.end(), edgeComparator);

    TourParts tourParts(tsplibProblem.getDimension());
    for (std::pair<unsigned int, unsigned int> edge : edges) {
        tourParts.join(edge.first, edge.second);
        if (tourParts.isTourClosable()) {
            return tourParts.closeTour();
        }
    }
}
