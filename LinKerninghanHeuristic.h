//
// Created by karl on 14.05.19.
//

#ifndef LINKERNINGHANALGORITHM_LINKERNINGHANHEURISTIC_H
#define LINKERNINGHANALGORITHM_LINKERNINGHANHEURISTIC_H

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

Tour linKerninghanHeuristic(const TsplibProblem &tsplibProblem, const Tour &startTour);

#endif //LINKERNINGHANALGORITHM_LINKERNINGHANHEURISTIC_H
