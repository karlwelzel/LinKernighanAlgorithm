//
// Created by Karl Welzel on 29/03/2019.
//

#ifndef LINKERNINGHANALGORITHM_TSPLIBUTILS_H
#define LINKERNINGHANALGORITHM_TSPLIBUTILS_H

std::string trim(const std::string &, const std::string &);

class TsplibProblem {
private:
    const char delimiter = ':';

    std::string name;
    std::string type;
    unsigned int dimension = 0;
    std::string edgeWeightType;
    std::string edgeWeightFormat;
    std::string nodeCoordType = "TWOD_COORDS";

    // if EDGE_WEIGHT_TYPE is EXPLICIT
    std::vector<std::vector<int>> matrix;  // The matrix described in the TSPLIB file

    // if EDGE_WEIGHT_TYPE is *_2D
    std::vector<std::vector<double>> coordinates; // Every entry is a 2d coordinate

    std::string interpretKeyword(const std::string &keyword, const std::string &value);

public:
    TsplibProblem();

    std::string readFile(std::ifstream &inputFile);


    const std::string &getName() const;

    const std::string &getType() const;

    unsigned int getDimension() const;


    int dist(unsigned int, unsigned int) const;
};

class TourVertex {
private:
    unsigned int NO_VERTEX = std::numeric_limits<unsigned int>::max();

    std::pair<unsigned int, unsigned int> neighbors;

public:
    TourVertex();

    TourVertex(unsigned int neighbor1, unsigned int neighbor2);

    explicit TourVertex(const std::pair<unsigned int, unsigned int> &neighbors);

    const std::pair<unsigned int, unsigned int> &getNeighbors() const;

    void setNeighbors(const std::pair<unsigned int, unsigned int> &neighbors);

    void setNeighbors(unsigned int neighbor1, unsigned int neighbor2);
};

// TODO: Refactor tour to be a map to individual TourVertex's that point to their neighbors
class Tour {
protected:
    // The vertices of this tour
    std::map<unsigned int, TourVertex> vertices;

    void setVertices(const std::list<unsigned int> &tourList);

public:
    Tour();

    explicit Tour(const std::list<unsigned int> &tourList);

    const TourVertex tourVertexAt(const unsigned int &i) const;

    const unsigned int next(unsigned int previous, unsigned int current) const;

    void setNext(unsigned int previous, unsigned int current, unsigned int next);

    const unsigned int length(TsplibProblem &tsplibProblem);
};

std::ostream &operator<<(std::ostream &, const Tour &);


class TsplibTour : public Tour {
private:
    const char delimiter = ':';

    std::string name;
    std::string type;
    unsigned int dimension = 0;

    std::string interpretKeyword(const std::string &keyword, const std::string &value);

public:
    const std::string &getName() const;

    const std::string &getType() const;

    std::string readFile(std::ifstream &inputFile);
};

#endif //LINKERNINGHANALGORITHM_TSPLIBUTILS_H
