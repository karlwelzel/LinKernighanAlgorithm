//
// Created by Karl Welzel on 14.05.19.
//

#ifndef LINKERNINGHANALGORITHM_LINKERNIGHANHEURISTIC_H
#define LINKERNINGHANALGORITHM_LINKERNIGHANHEURISTIC_H

#include <cstddef>
#include <functional>
#include <iostream>
#include <vector>
#include "Tour.h"
#include "TsplibUtils.h"

// ============================================ CandidateEdgeGraph class ===============================================

// This class provides candidate edges for each vertex in a graph and functions for generating them

class CandidateEdges {
private:
    std::vector<std::vector<vertex_t>> neighbors;

    // Find the k nearest neighbors in the set of vertices 0, ..., dimension by using the distCompare function that
    // decides for three vertices v, w1 and w2 if the distance between v and w1 is smaller than the distance between v
    // and w2
    static CandidateEdges rawNearestNeighbors(dimension_t dimension, std::size_t k,
                                              const std::function<bool(vertex_t, vertex_t, vertex_t)> &distCompare);

public:
    enum Type {
        ALL_NEIGHBORS, NEAREST_NEIGHBORS, ALPHA_NEAREST_NEIGHBORS, OPTIMIZED_ALPHA_NEAREST_NEIGHBORS
    };

    CandidateEdges() = default;

    // Fill edges with dimension copies of fillValue (just as the constructor of std::vector)
    CandidateEdges(dimension_t dimension, const std::vector<vertex_t> &fillValue);

    // For each vertex choose all edges as candidate edges
    static CandidateEdges allNeighbors(const TsplibProblem &problem);

    // For each vertex choose the k edges with minimal distance as candidate edges
    static CandidateEdges nearestNeighbors(const TsplibProblem &problem, std::size_t k);

    // For each vertex choose the k edges with minimal alpha distance as candidate edges
    // The alpha distance of an edge is defined as the increase in length of a 1-tree when required to contain this edge
    static CandidateEdges alphaNearestNeighbors(const TsplibProblem &problem, std::size_t k);

    // For each vertex choose the k edges with minimal alpha distance as candidate edges
    // The alpha distances are optimized with subgradient optimization, see AlphaDistance
    static CandidateEdges optimizedAlphaNearestNeighbors(const TsplibProblem &problem, std::size_t k);

    // Create candidate edges of type candidateEdgeType with k candidate edges for each vertex
    // k is ignored for Type::ALL_NEIGHBORS
    static CandidateEdges create(const TsplibProblem &problem, Type candidateEdgeType, std::size_t k);

    // Forwards the [] operator of neighbors
    std::vector<vertex_t> &operator[](std::size_t index);
};

// ========================================== LinKernighanHeuristic class ==============================================

// This class represents a single run of the Lin-Kernighan-heuristic. Each run consists of multiple trials and in every
// trial a random tour is improved until no further improvement is found. The generation of random tours and the search
// for improvements can depend on previous trials.

class LinKernighanHeuristic {
private:
    const std::size_t backtrackingDepth = 5;
    const std::size_t infeasibilityDepth = 2;

    // The TsplibProblem that should be solved
    TsplibProblem tsplibProblem;

    // The best solution tour found by the algorithm (initially none)
    Tour currentBestTour;

    // The candidate edges used
    CandidateEdges candidateEdges;

    // Chooses a random element from the vector elements
    static vertex_t chooseRandomElement(const std::vector<vertex_t> &elements);

    // Generates a random good tour based on the current best tour and the candidate edges
    Tour generateRandomTour();

    // The core part of the algorithm as described in Combinatorial Optimization
    Tour improveTour(const Tour &startTour);

public:
    LinKernighanHeuristic() = delete;

    explicit LinKernighanHeuristic(TsplibProblem &tsplibProblem, CandidateEdges candidateEdges);

    // Return the best tour found after numberOfTrials trials. If the relative increase of the length of the best tour
    // compared to optimumTourLength is below acceptableError the algorithm will stop and return it immediately.
    // verboseOutput turns debug output on or off
    Tour findBestTour(std::size_t numberOfTrials, distance_t optimumTourLength = 0, double acceptableError = 0,
                      bool verboseOutput = true);
};

#endif //LINKERNINGHANALGORITHM_LINKERNIGHANHEURISTIC_H
