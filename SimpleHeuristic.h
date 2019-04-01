//
// Created by Karl Welzel on 29/03/2019.
//

#ifndef LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H
#define LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H

class TourParts { // Enhanced UnionFind structure
private:
    std::vector<unsigned int> parent;
    std::vector<unsigned int> rank;
    std::vector<std::list<unsigned int>> tourPart;

public:
    explicit TourParts(unsigned int dimension);

    const std::list<unsigned int> &getTourPartAt(unsigned int i) const;

    unsigned int find(unsigned int x);

    int join(unsigned int x, unsigned int y);
};

Tour simpleHeuristic(TsplibProblem &tsplibProblem);

#endif //LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H
