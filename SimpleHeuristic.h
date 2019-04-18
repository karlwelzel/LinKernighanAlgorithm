//
// Created by Karl Welzel on 29/03/2019.
//

#ifndef LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H
#define LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H

class TourParts : public VertexList {
private:
    unsigned int dimension;
    std::vector<unsigned int> parent;
    std::vector<unsigned int> size;
    std::vector<std::pair<unsigned int, unsigned int>> pathEnds;

public:
    explicit TourParts(unsigned int dimension);

    unsigned int find(unsigned int x);

    void join(unsigned int x, unsigned int y);

    bool isTourClosable();

    Tour closeTour();
};

Tour simpleHeuristic(TsplibProblem &tsplibProblem);

#endif //LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H
