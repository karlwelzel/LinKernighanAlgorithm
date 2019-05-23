//
// Created by karl on 14.05.19.
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


Tour linKernighanHeuristic(const TsplibProblem &tsplibProblem, const Tour &startTour);

#endif //LINKERNINGHANALGORITHM_LINKERNIGHANHEURISTIC_H
