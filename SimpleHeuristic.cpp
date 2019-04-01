//
// Created by Karl Welzel on 29/03/2019.
//

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

TourParts::TourParts(unsigned int dimension) {
    for (unsigned int i = 0; i < dimension; ++i) {
        parent.push_back(i); // parent[i] = i
        rank.push_back(0);   // rank[i] = 0
        // The tour only consists of the vertex i
        tourPart.push_back(std::list<unsigned int>{i});
    }
}

const std::list<unsigned int> &TourParts::getTourPartAt(unsigned int i) const {
    return tourPart.at(i);
}

unsigned int TourParts::find(unsigned int x) {
    if (parent.at(x) != x) {
        parent.at(x) = find(parent.at(x));
    }
    return parent.at(x);
}

int TourParts::join(unsigned int x, unsigned int y) {
    // Join the parts of a tour that x and y are part of by adding the edge (x, y)
    unsigned int xr = find(x); // root of x
    unsigned int yr = find(y); // root of y
    if (xr == yr) {
        return -1; // A tour part cannot be joined with itself, this would lead to non-Hamiltonian tours
    }
    if (x == tourPart.at(xr).front() and y == tourPart.at(yr).front()) {
        tourPart.at(xr).reverse();
    } else if (x == tourPart.at(xr).front() and y == tourPart.at(yr).back()) {
        tourPart.at(xr).reverse();
        tourPart.at(yr).reverse();
    } else if (x == tourPart.at(xr).back() and y == tourPart.at(yr).front()) {
    } else if (x == tourPart.at(xr).back() and y == tourPart.at(yr).back()) {
        tourPart.at(yr).reverse();
    } else {
        return -1; // x and y cannot be joined, because one of them is not an end of its tour part
    }
    // The tour parts can now be joined by adding tourPart.at(yr) at the end of tourPart.at(xr) or by adding
    // tourPart.at(xr) at the beginning of tourPart.at(yr)
    unsigned int newRoot;
    if (rank.at(xr) > rank.at(yr)) {
        newRoot = xr;
        parent.at(yr) = xr;
        tourPart.at(xr).splice(tourPart.at(xr).end(), tourPart.at(yr));
    } else {
        newRoot = yr;
        parent.at(xr) = yr;
        tourPart.at(yr).splice(tourPart.at(yr).begin(), tourPart.at(xr));
    }
    if (rank.at(xr) == rank.at(yr)) {
        rank.at(yr)++;
    }
    return newRoot;
}

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


Tour simpleHeuristic(TsplibProblem &tsplibProblem) {
    std::vector<std::pair<unsigned int, unsigned int>> edges;
    for (int i = 0; i < tsplibProblem.getDimension(); ++i) {
        for (int j = i + 1; j < tsplibProblem.getDimension(); ++j) {
            edges.emplace_back(i, j);
        }
    }

    EdgeCostComparator edgeComparator(tsplibProblem);
    std::sort(edges.begin(), edges.end(), edgeComparator);

    TourParts tourParts(tsplibProblem.getDimension());
    for (std::pair<unsigned int, unsigned int> edge : edges) {
        int root = tourParts.join(edge.first, edge.second);
        if (root >= 0 and tourParts.getTourPartAt(root).size() == tsplibProblem.getDimension()) {
            return Tour(tourParts.getTourPartAt(root));
        }
    }
}
