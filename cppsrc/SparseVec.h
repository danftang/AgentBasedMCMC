//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_SPARSEVEC_H
#define GLPKTEST_SPARSEVEC_H


#include <cstdlib>
#include <ostream>
#include <misc/fvs.h>

class SparseVec : public FVS {
public:

    SparseVec(int dimension);
    ~SparseVec();

    // zero-index based accessors
    double operator [](int i) const;
    double &operator [](int i);

    void toDense(double *denseVec) const;
    void entry(int k, std::pair<int,double> &);
    void add(int i, double v);
    void clear();
    double dotProd(double *dense) const;

protected:
};

std::ostream &operator <<(std::ostream &out, const SparseVec &sVector);

#endif //GLPKTEST_SPARSEVEC_H
