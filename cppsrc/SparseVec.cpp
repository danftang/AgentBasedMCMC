//
// Created by daniel on 09/04/2021.
//

#include <iomanip>
#include "SparseVec.h"

SparseVec::SparseVec(int dim) {
    sparseSize = 0;
    dimension = dim;
    indices = new int[dimension+1];
    vals = new double[dimension+1];
}

SparseVec::~SparseVec() {
    free(indices);
    free(vals);
}

double &SparseVec::operator[](int i) {
    for(int k=1; k<=sparseSize; k++) {
        if(indices[k] == i) return vals[k];
    }
    sparseSize++;
    indices[sparseSize] = i;
    return vals[sparseSize]; // allows zero values
}


double SparseVec::operator[](int i) const {
    for(int k=1; k<=sparseSize; k++) {
        if(indices[k] == i) return vals[k];
    }
    return 0.0;
}

void SparseVec::add(int i, double v) {
    indices[++sparseSize] = i;
    vals[sparseSize] = v; // allows zero values
}

std::ostream &operator <<(std::ostream &out, SparseVec &sVector) {
    int i;
    double dense[sVector.dimension];

    for(i=0; i<sVector.dimension; ++i) { dense[i] = 0.0;}
    for(i=1; i<=sVector.sparseSize; ++i) {
        if(sVector.indices[i]-1 >= sVector.dimension) out << "Out of range index!!";
        dense[sVector.indices[i]-1] = sVector.vals[i];
    }

    for(i=0; i<sVector.dimension; ++i) {
        out << std::setw(12) << dense[i] << "\t";
    }
    return out;
}
