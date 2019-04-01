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


class Tour {
protected:
    std::list<unsigned int> vertices; // The vertices of this tour in the correct order

public:
    Tour();

    explicit Tour(const std::list<unsigned int> &tour);


    const std::list<unsigned int> &getVertices() const;


    const std::list<unsigned int>::iterator iteratorOf(unsigned int vertex);

    const std::list<unsigned int>::iterator after(std::list<unsigned int>::iterator it);

    const std::list<unsigned int>::iterator before(std::list<unsigned int>::iterator it);

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
