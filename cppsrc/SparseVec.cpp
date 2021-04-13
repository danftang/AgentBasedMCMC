//
// Created by daniel on 09/04/2021.
//

#include <iomanip>
#include <iostream>
#include "SparseVec.h"

SparseVec::SparseVec(int dim) {
    n = 0;
    nnz = dim;
    ind = new int[dim] -1;
    vec = new double[dim] -1;
}

SparseVec::~SparseVec() {
    free(ind+1);
    free(vec+1);
}

double &SparseVec::operator[](int i) {
    for(int k=1; k<=nnz; k++) {
        if(ind[k] == i+1) return vec[k];
    }
    ind[++nnz] = i+1;
    return vec[nnz]; // allows zero values
}


double SparseVec::operator[](int i) const {
    for(int k=1; k<nnz; k++) {
        if(ind[k] == i+1) return vec[k];
    }
    return 0.0;
}

void SparseVec::add(int i, double v) {
    if(v != 0.0) {
        ind[++nnz] = i + 1;
        vec[nnz] = v;
    } // doesn't delete if index already exists
}

void SparseVec::clear() {
    nnz = 0;
}

void SparseVec::toDense(double *dense) const {
    int i;
    for(i=0; i<n; ++i) { dense[i] = 0.0; }
    for(i=1; i<=nnz; ++i) {
        if(ind[i] > n || ind[i] < 1)
            std::cout << "Out of range index[" << i << "] = " << ind[i] << " -> " << vec[i] << std::endl;
        dense[ind[i]-1] = vec[i];
    }
}

void SparseVec::entry(int k, std::pair<int, double> &retEntry) {
    int k1 = k+1;
    retEntry.first = ind[k1];
    retEntry.second = vec[k1];
}

double SparseVec::dotProd(double *dense) const {
    double *dense1base = dense-1;
    double dp = 0.0;
    int j;
    for(int i=1; i<=nnz; ++i) {
        j = ind[i];
        dp += dense1base[j] * vec[i];
    }
    return dp;
}


std::ostream &operator <<(std::ostream &out, const SparseVec &sVector) {
    int i;
    double dense[sVector.n];
    sVector.toDense(dense);
    for(i=0; i<sVector.n; ++i) out << std::setw(12) << dense[i] << "\t";
    return out;
}
