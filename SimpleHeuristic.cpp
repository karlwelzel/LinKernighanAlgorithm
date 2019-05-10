//
// Created by Karl Welzel on 29/03/2019.
//

// This library contains the code for the introduction assignment

#include <utility>
#include <tuple>
#include <vector>
#include <algorithm>
#include "SimpleHeuristic.h"


// ============================================== TourParts class ======================================================

TourParts::TourParts(const dimension_t dimension) : VertexList(dimension) {
    for (dimension_t i = 0; i < dimension; ++i) {
        parent.push_back(i); // parent[i] = i
        size.push_back(1);   // size[i] = 1
        pathEnds.emplace_back(i, i);
    }
}

vertex_t TourParts::find(const vertex_t x) {
    if (parent.at(x) != x) { // path compression
        parent.at(x) = find(parent.at(x));
    }
    return parent.at(x);
}

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

bool TourParts::isTourClosable() {
    return size.at(find(0)) == getDimension();
}

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


Tour simpleHeuristic(const TsplibProblem &tsplibProblem) {
    std::vector<std::tuple<distance_t, vertex_t, vertex_t>> edges;
    for (vertex_t i = 0; i < tsplibProblem.getDimension(); ++i) {
        for (vertex_t j = i + 1; j < tsplibProblem.getDimension(); ++j) {
            edges.emplace_back(tsplibProblem.dist(i, j), i, j);
        }
    }

    // Sort edges by edge cost/distance in ascending order
    std::sort(edges.begin(), edges.end());

    TourParts tourParts(tsplibProblem.getDimension());
    for (const auto edge : edges) {
        vertex_t vertex1, vertex2;
        std::tie(std::ignore, vertex1, vertex2) = edge;
        tourParts.join(vertex1, vertex2);
        if (tourParts.isTourClosable()) {
            return tourParts.closeTour();
        }
    }

    // This line will never be reached, because the graph is a complete graph
    throw std::runtime_error("SimpleHeuristic: Although every edge was checked, no tour has been found.");
}
