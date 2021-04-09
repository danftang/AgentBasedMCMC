//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_SPARSEVEC_H
#define GLPKTEST_SPARSEVEC_H


#include <cstdlib>
#include <ostream>

class SparseVec {
public:
    int         dimension;
    int         sparseSize;
    int *       indices;
    double *    vals;

    SparseVec(int dimension);
    ~SparseVec();

    double operator [](int i) const;
    double &operator [](int i);

    void add(int i, double v);
};

std::ostream &operator <<(std::ostream &out, SparseVec &sVector);

#endif //GLPKTEST_SPARSEVEC_H
