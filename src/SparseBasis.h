// Represents a lattice of points defined by a basis in an "augmented" space, stored in sparse format
// The valid points of the lattice are contained within a hyperrectangle with one corner at the origin
// and an diagonally opposite corner at "H", so that the valid points can be defined as the set
//
// { X : X = BX' - F, 0 <= X <= H }
//
// where the columns of B are the basis vectors, and the elements of $X'$ are integers.
//
// Created by daniel on 07/03/2022.
//

#ifndef ABMCMC_SPARSEBASIS_H
#define ABMCMC_SPARSEBASIS_H

#include <vector>
#include "TableauNormMinimiser.h"


template<class T>
class SparseBasis: public std::vector<SparseVec<T>> {
public:
    std::vector<double>     H; // The upper bounds of the hyperrectangle (see intro notes above)
    std::vector<double>     X; // the current point on the lattice.

    SparseBasis(TableauNormMinimiser &tableau);
    SparseBasis(const ConvexPolyhedron &polyhedron);

    int nRows() const { return H.size(); }
    int nCols() const { return this->size(); }
//    friend std::ostream &operator <<(std::ostream &out, const SparseBasis &basis);
};


template<class T>
SparseBasis<T>::SparseBasis(TableauNormMinimiser &tableau):
        std::vector<SparseVec<T>>(tableau.nNonBasic()),
        H(tableau.cols.size()+tableau.nAuxiliaryVars),
        X(H.size())
{
    int jBasis = 0;
    for(int j=0; j<tableau.cols.size(); ++j) {
        if(!tableau.cols[j].isBasic) {
            (*this)[jBasis].reserve(tableau.cols[j].size()+1);
            for(auto i : tableau.cols[j]) {
                int basisi = tableau.minimalBasis[i];
                if(basisi < 0) basisi = -basisi + tableau.cols.size(); // transform auxiliaries to after end of X
                (*this)[jBasis].insert(basisi, tableau.rows[i].at(j));
            }
            (*this)[jBasis].insert(j,1); // add the non-basic itself
            ++jBasis;
        }
    }

}

template<class T>
SparseBasis<T>::SparseBasis(const ConvexPolyhedron &polyhedron) {
    TableauNormMinimiser tableau(polyhedron);

    H.resize(tableau.cols.size()+tableau.nAuxiliaryVars);
    X.resize(H.size());
    this->resize(tableau.nNonBasic());

    // transfer basis over from tableau
    int jBasis = 0;
    for(int j=0; j<tableau.cols.size(); ++j) {
        if(!tableau.cols[j].isBasic) {
            (*this)[jBasis].reserve(tableau.cols[j].size()+1);
            for(int i : tableau.cols[j]) {
                int basisi = tableau.minimalBasis[i];
                if(basisi < 0) basisi = -basisi + tableau.cols.size(); // transform auxiliaries to after end of X
                (*this)[jBasis].insert(basisi, tableau.rows[i].at(j));
            }
            (*this)[jBasis].insert(j,1); // add the non-basic itself
            ++jBasis;
        }
    }

    // initialise H and X
    for(int i=0; i<H.size(); ++i) {

    }

}

template<class T>
std::ostream &operator <<(std::ostream &out, const SparseBasis<T> &basis) {
    for(int i=0; i<basis.nRows(); ++i) {
        for(int j=0; j<basis.nCols(); ++j) {
            out << basis[j][i] << "\t";
        }
        out << std::endl;
    }
    return out;
}



#endif //ABMCMC_SPARSEBASIS_H
