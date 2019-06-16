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
    enum Type {
        ALL_NEIGHBORS, NEAREST_NEIGHBORS, ALPHA_NEAREST_NEIGHBORS
    };

    // For each vertex choose all edges as candidate edges
    static CandidateEdges allNeighbors(const TsplibProblem &problem);

    // For each vertex choose the k edges with minimal distance as candidate edges
    static CandidateEdges nearestNeighbors(const TsplibProblem &problem, size_t k = 10);

    // For each vertex choose the k edges with minimal alpha distance as candidate edges
    // The alpha distance of an edge is defined as the increase in length of a 1-tree when required to contain this edge
    static CandidateEdges alphaNearestNeighbors(const TsplibProblem &problem, size_t k = 10);
};

// ========================================== LinKernighanHeuristic class ==============================================

// This class represents a single run of the Lin-Kernighan-heuristic. Each run consists of multiple trials and in every
// trial a random tour is improved until no further improvement is found. The generation of random tours and the search
// for improvements can depend on previous trials.

class LinKernighanHeuristic {
private:
    const size_t backtrackingDepth = 5;
    const size_t infeasibilityDepth = 2;

    // The TsplibProblem that should be solved
    TsplibProblem tsplibProblem;

    // The best solution tour found by the algorithm (initially none)
    Tour currentBestTour;

    // The candidate edges used
    CandidateEdges candidateEdges;

    static vertex_t chooseRandomElement(const std::vector<vertex_t> &elements);

    // Generates a random good tour based on the current best tour and the candidate edges
    Tour generateRandomTour();

    // The core part of the algorithm as described in Combinatorial Optimization
    Tour improveTour(const Tour &startTour);

public:
    LinKernighanHeuristic() = delete;

    explicit LinKernighanHeuristic(TsplibProblem &tsplibProblem,
                                   CandidateEdges::Type candidateEdgeType = CandidateEdges::ALPHA_NEAREST_NEIGHBORS);

    Tour findBestTour(size_t numberOfTrials = 10);
};

#endif //LINKERNINGHANALGORITHM_LINKERNIGHANHEURISTIC_H
