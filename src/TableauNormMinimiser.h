//
// Created by daniel on 10/12/2021.
//

#ifndef GLPKTEST_TABLEAUNORMMINIMISER_H
#define GLPKTEST_TABLEAUNORMMINIMISER_H

#include <set>
#include <vector>
#include <map>
#include <list>
#include "ConvexPolyhedron.h"

// Takes an LP Problem and finds the basis that pivots out as many
// fixed-variables as possible, while attempting to keep the L0-norm of the tableau
// as low as possible.
//
// Starting with a set of linear constraints on X
// L <= CX <= U, DX = E, 0 <= X <= H
// we introduce auxiliary variables, A, so that
// (C -I)(X A)^T = (L E)^T,  0 <= (X A) <= (H U-L)
// (D  0)
// If D has M rows, we can choose any M linearly independent columns
// and pivot so that they have only one non-zero element equal to -1.
// to give
// M(X A)^T + F = 0,  0 <= (X A) <= (H U-L)
// where M has a total of |A| + M reduced columns (i.e. basic variables)
//
// Given this, we can transform into the form
// Q X_N + F = (X A), 0 <= (X A) <= (H U-L)
// by taking basic variables over to the right-hand side and adding rows to Q
// to insert the non-basic vars.
// The columns of Q form the sparse basis so if M is sparse, so is Q.
//
// All fixed variables are assumed to be functions of integers only.
// Pivoting only occurs on elements that have value +-1 to retain
// integer coefficients.
//
// TODO: Need to maintain lower limits so that factors get correct values after
//  reducing the space
//
template<class T>
class TableauNormMinimiser {
public:

    static constexpr int maxVectorsToTest = 10;


    class Column: public std::set<int> { // set of non-zero orw indices
    public:
        std::list<int>::iterator sparsityEntry; // iterator into colsBySaprsity entry
        bool isBasic;

//        explicit Column(const SparseVec<double> &col);
        Column(): isBasic(false) {};
    };

    class Row: public std::map<int,T> { // map from non-zero col index to element value
    public:
        std::list<int>::iterator sparsityEntry; // iterator into rowsBySaprsity entry
        bool isActive;

        Row(const SparseVec<T> &row, bool isActive);

        Row &operator *=(T c);
    };

    std::vector<Column>     cols;   // tableau columns
    std::vector<Row>        rows;   // tableau rows

    std::vector<int>            constraintIndexByCol;   // index of rectilinear constraints in original problem
    std::vector<int>            constraintIndexByAuxiliary;   // index of auxiliary constraint in the original problem
    std::vector<std::list<int>> colsBySparsity;     // only contains non-basic cols
    std::vector<std::list<int>> rowsBySparsity;     // only contains active rows
    std::vector<int>            basis;            // basic variables: -ve means auxiliary, otherwise col index
    int                         nAuxiliaryVars;     // number of auxiliary variables
    std::vector<T>              F;                  // aggregated value of all non-basic fixed vars
    std::vector<T>              Hc;                  // upper bounds of columns
    std::vector<T>              Ha;                  // upper bounds of auxiliaries

//    TableauNormMinimiser(glp::Problem &problem);
    TableauNormMinimiser(const ConvexPolyhedron<T> &problem);

    int nNonBasic() const { return cols.size() + nAuxiliaryVars - rows.size(); }

    void findMinimalBasis();

    double operator ()(int i, int j) const { auto it = rows[i].find(j); return it == rows[i].end()?0.0:it->second; }

    std::pair<int,int> findMarkowitzPivot();

    void pivot(int i,int j);

    int sparsestColInRow(int i);
    int sparsestRowInCol(int j);

//    void updateRowSparsity(int i);
    void addRowSparsityEntry(int i);
    void removeRowSparsityEntry(int i);
//    void updateColSparsity(int j);
    void addColSparsityEntry(int j);
    void removeColSparsityEntry(int j);

    void setColBasic(int j, int i);
    void inactivateRow(int i);

    double meanColumnL0Norm() const;
    double meanColumnL1Norm() const;

};

template<class T>
std::ostream &operator <<(std::ostream &out, const TableauNormMinimiser<T> &tableau) {
    // show upper limits
    out << "    \t \t";
    for(int j=0; j<tableau.cols.size(); ++j) {
        out << tableau.Hc[j] << "\t";
    }
    out << std::endl;
    // print tableau
    for(int i=0; i<tableau.rows.size(); ++i) {
        out << "x(" << tableau.basis[i] << ")\t=\t";
        for(int j=0; j<tableau.cols.size(); ++j) {
            out << tableau(i,j) << "\t";
        }
        out << "+\t" << tableau.F[i];
        if(tableau.basis[i] < 0) out << " <= " << tableau.Ha[-tableau.basis[i]];
        out << std::endl;
    }
    // mark basic columns
    out << "    \t \t";
    for(int j=0; j<tableau.cols.size(); ++j) {
        out << (tableau.cols[j].isBasic?"B\t":"-\t");
    }
    out << std::endl;
    return out;
}

// Turns the supplied convex polyhedron into a set of equality constraints by introducing auxiliary
// variables where the upper and lower bounds are not equal. So if we have
// L <= CX <= U
// then
// (-L M)(1 X A)^T = 0, 0 <= A <= U-L
// where the zero'th elelemet of (1 X A) is 1.
template<class T>
TableauNormMinimiser<T>::TableauNormMinimiser(const ConvexPolyhedron<T> &problem)
{
    basis.reserve(problem.size());
    F.reserve(problem.size());
    rows.reserve(problem.size());
    Ha.reserve(problem.size());
    Ha.push_back(0.0);
    constraintIndexByAuxiliary.reserve(problem.size());
    constraintIndexByAuxiliary.push_back(-1);
    int maxCol = 0; // maximum column id seen so far
    nAuxiliaryVars = 0;
    int constraintIndex = 0;
    for(const Constraint<T> &constraint: problem) {
        if(constraint.coefficients.sparseSize() == 1) { // rectilinear limits, don't add row
            assert(constraint.lowerBound == 0);
            int xi = constraint.coefficients.indices[0];
            if(Hc.size() < xi + 1) Hc.resize(xi+1);
            Hc[xi] = constraint.upperBound/constraint.coefficients.values[0];
            if(constraintIndexByCol.size() < xi+1) constraintIndexByCol.resize(xi+1,-1);
            constraintIndexByCol[xi] = constraintIndex;
            if(xi > maxCol) maxCol = xi;
        } else { // normal tableau entry
            bool isActive = (constraint.upperBound == constraint.lowerBound);
            rows.push_back(Row(constraint.coefficients, isActive));
            if (isActive) {
                addRowSparsityEntry(rows.size()-1);
                basis.push_back(0); // fixed var
            } else {
                basis.push_back(-(++nAuxiliaryVars));
                Ha.push_back(constraint.upperBound - constraint.lowerBound);
                constraintIndexByAuxiliary.push_back(constraintIndex);
            }
            F.push_back(-constraint.lowerBound); // lower bounds are all set to 0
            int rowMaxCol = rows.back().rbegin()->first;
            if(rowMaxCol > maxCol) maxCol = rowMaxCol;
        }
        ++constraintIndex;
    }

    Hc.resize(maxCol+1);
    cols.resize(maxCol +1);
    constraintIndexByCol.resize(maxCol + 1,-1);
    for(int j=1; j<=maxCol; ++j) assert(Hc[j] != 0); // no rectilinear limit given for non-basic j or fixed at 0
    for(int i=0; i<rows.size(); ++i) {
        for(auto [j, v]: rows[i]) {
            cols[j].insert(i);
        }
    }
    for(int j=0; j<=maxCol; ++j) {
        addColSparsityEntry(j);
    }
    std::cout << "Constructed tableau:" << std::endl << *this << std::endl;
    std::cout << "Initial mean L0 norm = " << meanColumnL0Norm() << std::endl;
    std::cout << "Initial L1 norm = " << meanColumnL1Norm() << std::endl;
    findMinimalBasis();
}

template<class T>
void TableauNormMinimiser<T>::findMinimalBasis() {
    while(!rowsBySparsity.empty()) {
        auto [i,j] = findMarkowitzPivot();
        pivot(i,j);
//        std::cout << minimalBasis << std::endl;
//        std::cout << rowsBySparsity << std::endl;
    }
    std::cout << "Minimised tableau:" << std::endl << *this << std::endl;
    std::cout << "Reduced mean L0 norm = " << meanColumnL0Norm() << std::endl;
    std::cout << "Reduced L1 norm = " << meanColumnL1Norm() << std::endl;

}

template<class T>
std::pair<int,int> TableauNormMinimiser<T>::findMarkowitzPivot() {
    int k=1;
    int lowestScore = INT32_MAX;
    int vectorsTested = 0;
    std::pair<int,int> bestPivot(-1,-1);
    while(k<colsBySparsity.size() && k<rowsBySparsity.size()) {
//        std::cout << "k = " << k << std::endl;
        for(int j: colsBySparsity[k]) {
            int bestRow = sparsestRowInCol(j);
            if(bestRow != -1) {
                int score = (rows[bestRow].size() - 1) * (cols[j].size() - 1);
                if (score < lowestScore) {
                    bestPivot = std::pair<int, int>(bestRow, j);
                    if (score <= (k - 1) * (k - 1)) return bestPivot;
                    lowestScore = score;
                }
                if (++vectorsTested >= maxVectorsToTest) return bestPivot;
            }
        }
        for(int i: rowsBySparsity[k]) {
            int bestCol = sparsestColInRow(i);
            if(bestCol != -1) {
                int score = (rows[i].size() - 1) * (cols[bestCol].size() - 1);
                if (score < lowestScore) {
                    bestPivot = std::pair<int, int>(i, bestCol);
                    if (score <= (k - 1) * (k - 1)) return bestPivot;
                    lowestScore = score;
                }
                if (++vectorsTested >= maxVectorsToTest) return bestPivot;
            }
        }
        ++k;
    }
    return bestPivot;
}


template<class T>
void TableauNormMinimiser<T>::pivot(int pi, int pj) {
    assert(rows[pi].isActive);
//    std::cout << "Pivoting at " << pi << " " << pj << std::endl;
    T Mpipj = rows[pi].at(pj);

    // remove sparsity entries of all affected rows and columns
    for( auto [j, v] : rows[pi]) removeColSparsityEntry(j);

    rows[pi] *= -1.0/Mpipj;     // divide row through to make the pivot point -1
    F[pi] *= -1.0/Mpipj;

    for(int i : cols[pj]) {
        if(i != pi) {
            if(rows[i].isActive) removeRowSparsityEntry(i);
            T Mipj = rows[i].at(pj);
            for (const auto[j, Mpij]: rows[pi]) {
                if(j != pj) {
                    auto MijIt = rows[i].find(j);
                    if (MijIt == rows[i].end()) { // fill-in
                        rows[i][j] = Mipj * Mpij;
                        cols[j].insert(i);
                    } else {
                        T newMij = (rows[i][j] += Mipj * Mpij);
                        if (newMij == 0.0) { // drop-out
                            cols[j].erase(i);
                            rows[i].erase(j);
                        }
                    }
                }
            }
            F[i] += Mipj * F[pi];
        }
    }

    for(int i : cols[pj]) {
        if(i != pi) {
            rows[i].erase(pj);
            if(rows[i].isActive) addRowSparsityEntry(i);
        }
    }
    setColBasic(pj, pi);
    for(const auto [j, Mpij]: rows[pi]) {
        if(j != pj) addColSparsityEntry(j);
    }
    inactivateRow(pi);
}


//void TableauNormMinimiser::updateRowSparsity(int i) {
//    Row &row = rows[i];
//    if(row.sparseSize != row.size()) {
//        removeRowSparsityEntry(i);
//        addRowSparsityEntry(i);
//    }
//}



template<class T>
void TableauNormMinimiser<T>::addRowSparsityEntry(int i) {
    Row &row = rows[i];
    if(rowsBySparsity.size() <= row.size()) rowsBySparsity.resize(row.size()+1);
    rowsBySparsity[row.size()].push_front(i);
    row.sparsityEntry = rowsBySparsity[row.size()].begin();
}

template<class T>
void TableauNormMinimiser<T>::removeRowSparsityEntry(int i) {
    Row &row = rows[i];
    assert(row.isActive);
    rowsBySparsity[row.size()].erase(row.sparsityEntry);
    while(rowsBySparsity.size() > 0 && rowsBySparsity.back().empty()) rowsBySparsity.pop_back();
}

template<class T>
void TableauNormMinimiser<T>::addColSparsityEntry(int j) {
    Column &col = cols[j];
//    col.sparseSize = col.size();
    if(colsBySparsity.size() <= col.size()) colsBySparsity.resize(col.size()+1);
    colsBySparsity[col.size()].push_front(j);
    col.sparsityEntry = colsBySparsity[col.size()].begin();
}

template<class T>
void TableauNormMinimiser<T>::removeColSparsityEntry(int j) {
    Column &col = cols[j];
    assert(!col.isBasic);
    colsBySparsity[col.size()].erase(col.sparsityEntry);
    col.sparsityEntry = colsBySparsity[col.size()].end();
    while(colsBySparsity.size() > 0 && colsBySparsity.back().empty()) colsBySparsity.pop_back();
}


template<class T>
TableauNormMinimiser<T>::Row::Row(const SparseVec<T> &row, bool isActive): isActive(isActive) {
    for (int k = 0; k < row.sparseSize(); ++k) {
        (*this)[row.indices[k]] = row.values[k];
    }
}

template<class T>
typename TableauNormMinimiser<T>::Row &TableauNormMinimiser<T>::Row::operator*=(T c) {
    for(auto it = this->begin(); it != this->end(); ++it) it->second *= c;
    return *this;
}


//TableauNormMinimiser::Column::Column(const SparseVec<double> &col): isBasic(false) {
//    for (int k = 0; k < col.sparseSize(); ++k) {
//        insert(col.indices[k]);
//    }
//}


template<class T>
void TableauNormMinimiser<T>::setColBasic(int j, int i) {
    Column &col = cols[j];
    col.isBasic = true;
    col.clear();
    col.insert(i);
    basis[i] = j;
}

template<class T>
void TableauNormMinimiser<T>::inactivateRow(int i) {
    Row &row = rows[i];
    removeRowSparsityEntry(i);
    row.isActive = false;
}

template<class T>
int TableauNormMinimiser<T>::sparsestColInRow(int i) {
    int minColSparsity = INT32_MAX;
    int minj = -1;
    for(const auto [j, Mij] : rows[i]) {
        int colSize = cols[j].size();
        if(colSize < minColSparsity && !cols[j].isBasic && fabs(Mij) == 1.0) {
            minColSparsity = colSize;
            minj = j;
        }
    }
//    std::cout << "Considering pivot from row " << i << " " << minj << std::endl;
    return minj;
}

template<class T>
int TableauNormMinimiser<T>::sparsestRowInCol(int j) {
    int minRowSparsity = INT32_MAX;
    int mini = -1;
    for(int i : cols[j]) {
        int rowSize = rows[i].size();
        if(rowSize < minRowSparsity && rows[i].isActive && fabs(rows[i].at(j)) == 1.0) {
            minRowSparsity = rowSize;
            mini = i;
        }
    }
//    std::cout << "Considering pivot " << mini << " from col " << j << std::endl;
    return mini;
}

template<class T>
double TableauNormMinimiser<T>::meanColumnL0Norm() const {
    double normSum = 0.0;
    int nNonBasic = 0;
    for(int j=0; j<cols.size(); ++j) {
        if(!cols[j].isBasic) {
            normSum += cols[j].size();
            ++nNonBasic;
        }
    }
    return normSum / nNonBasic;
}

template<class T>
double TableauNormMinimiser<T>::meanColumnL1Norm() const {
    double normSum = 0.0;
    int nNonBasic = 0;
    for(int j=0; j<cols.size(); ++j) {
        if(!cols[j].isBasic) ++nNonBasic;
    }
    for(int i=0; i<rows.size(); ++i) {
        for(const auto [j, Mij] : rows[i]) {
            if(!cols[j].isBasic) normSum += fabs(Mij);
        }
    }
    return normSum / nNonBasic;
}


#endif //GLPKTEST_TABLEAUNORMMINIMISER_H
