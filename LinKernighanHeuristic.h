//
// Created by Karl Welzel on 14.05.19.
//

#ifndef LINKERNINGHANALGORITHM_LINKERNIGHANHEURISTIC_H
#define LINKERNINGHANALGORITHM_LINKERNIGHANHEURISTIC_H

#include <vector>
#include "Tour.h"
#include "TsplibUtils.h"

// ============================================= AlternatingWalk class =================================================

// This class represents an alternating walk and adds some utility functions to the vector class

class AlternatingWalk : public std::vector<vertex_t> {
public:
    // Adds the first vertex to the end, to get a closed alternating walk, and returns the result.
    // Does not change the walk it is called upon
    AlternatingWalk close() const;

    // Appends vertex at the end, closes the walk (see close) and returns the result
    // Does not change the walk it is called upon
    AlternatingWalk appendAndClose(vertex_t vertex) const;

    bool containsEdge(vertex_t vertex1, vertex_t vertex2) const;
};

// Output the walk to the stream out, typically used with std::cout
std::ostream &operator<<(std::ostream &out, const AlternatingWalk &walk);


// ============================================ CandidateEdgeGraph class ===============================================

// This class provides candidate edges for each vertex in a graph and functions for generating them

class CandidateEdges : public std::vector<std::vector<vertex_t>> {
public:
    // For each vertex choose all edges as candidate edges
    static CandidateEdges allNeighbors(const TsplibProblem &problem);

    // For each vertex choose the k edges with minimal distance as candidate edges
    static CandidateEdges nearestNeighbors(const TsplibProblem &problem, size_t k = 5);

    // For each vertex choose the k edges with minimal alpha distance as candidate edges
    // The alpha distance of an edge is defined as the increase in length of a 1-tree when required to contain this edge
    static CandidateEdges alphaNearestNeighbors(const TsplibProblem &problem, size_t k = 5);
};


// ============================================= linKernighanHeuristic =================================================

// The algorithm is implemented as described in Combinatorial Optimization

Tour linKernighanHeuristic(const TsplibProblem &tsplibProblem, const Tour &startTour, CandidateEdges &candidateEdges,
                           size_t backtrackingDepth = 5, size_t infeasibilityDepth = 2);

#endif //LINKERNINGHANALGORITHM_LINKERNIGHANHEURISTIC_H
