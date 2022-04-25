//
// Created by daniel on 10/12/2021.
//

#ifndef GLPKTEST_TABLEAUNORMMINIMISER_H
#define GLPKTEST_TABLEAUNORMMINIMISER_H

#include <set>
#include <vector>
#include <map>
#include <list>
#include "boost/serialization/map.hpp"
#include "boost/serialization/set.hpp"
#include "ConvexPolyhedron.h"

// Takes an LP Problem and finds the basis that pivots out as many
// fixed-variables as possible, while attempting to keep the L0-norm of the tableau
// as low as possible.
//
// Starting with a set of linear constraints on X
// L <= CX <= U, DX = E, G <= X <= H
// we introduce auxiliary variables, A, so that
// (C)(X) + ( 0) = (A), L <= A <= U, G <= X <= H
// (D)      (-E)   (0)
//
// If D has m rows, this class chooses m linearly independent columns
// and pivots on rows in D so that M is left with all zeroes in those columns
// and D is left as a (truncated) triangular with "diagonal" elements equal to -1
// to give
//
// (M)(X) + F = (A), L <= A <= U, G <= X <= H
// (J)          (0)
//
// where M has a total of m reduced columns. The reduced columns make the
// corresponding elements of X into basic variables, the remainder being
// non-basic.
//
// Given an assignment of values to non-basic variables, X_N, the basic variables,
// X_B, and the auxilliary variables, A, are uniquely defined and can be found by
// rearranging the columns of M and J so that
//
// (Q  0)(X_N) + F = (A), L <= A <= U, G_N <= X_N <= H_N, G_B <= X_B <= H_B
// (R  T)(X_B)       (0)
//
// So that
//
// (Q)X_N + F = (A)   , L <= A <= U, G_N <= X_N < H_N, G_B <= X_B <= H_B
// (R)          (X_B)
//
// where Q is M with the zero columns removed, X_N is X with the corresponding
// basic variables removed, R is J with basic cols removed and X_B are the
// basic variables.
//
// All fixed variables are assumed to be functions of integers only.
// Pivoting only occurs on elements that have value +-1 to retain
// integer coefficients.
//
// The non-basic columns of M and J form the sparse basis, so we want to pivot in
// such a way as to maintain the sparsity of these matrices. We do this using a
// greedy algorithm based on the Markovitz criterion. See
//  Maros, I., 2002, "Computational techniques of the sipmlex method", Springer
//
// TODO: Need to maintain lower limits so that factors get correct values after
//  reducing the space
//
template<class T>
class TableauNormMinimiser {
public:

    static constexpr int maxVectorsToTest = 10;


    class Column: public std::set<int> { // set of non-zero row indices (only non-reduced rows are stored)
    public:
        std::list<int>::iterator sparsityEntry; // iterator into colsBySaprsity entry
        bool isBasic;

//        explicit Column(const SparseVec<double> &col);
        Column(): isBasic(false) {};
    private:
        friend class boost::serialization::access;

        template <typename Archive>
        void serialize(Archive &ar, const unsigned int version) {
            ar & static_cast<std::set<int> &>(*this) & isBasic;
            ar & isBasic;
        }
    };

    class Row: public std::map<int,T> { // map from non-zero col index to element value
    public:
        std::list<int>::iterator sparsityEntry; // iterator into rowsBySaprsity entry
        bool isActive;                          // true if this row is an unreduced equality constraint

        Row(const SparseVec<T> &row, bool isActive);
        Row() {}

        Row &operator *=(T c);
    private:
        friend class boost::serialization::access;

        template <typename Archive>
        void serialize(Archive &ar, const unsigned int version) {
            ar & static_cast<std::map<int,T> &>(*this) & isActive;
        }
    };

    std::vector<Column>     cols;   // tableau columns
    std::vector<Row>        rows;   // tableau rows

    std::vector<std::list<int>> colsBySparsity;     // cols sorted by sparsity, only contains non-basic cols
    std::vector<std::list<int>> rowsBySparsity;     // rows sorted by sparsity, only contains active rows
    std::vector<int>            basis;            // basic variables by row. -ve means auxiliary, otherwise col index of basic var
    int                         nAuxiliaryVars;     // number of auxiliary variables
    std::vector<T>              F;                  // constant in linear equation (see intro above)
    std::vector<T>              L;                  // lower bounds by row (equalities are converted to zero)
    std::vector<T>              U;                  // upper bounds by row (equalities are converted to zero)
//    TableauNormMinimiser(glp::Problem &problem);
    TableauNormMinimiser(const ConvexPolyhedron<T> &problem);

    TableauNormMinimiser()=default;

    int nNonBasic() const { return cols.size() + nAuxiliaryVars - rows.size(); }

    bool isEqualityConstraint(int i) const { return L[i] == U[i]; }

    void findMinimalBasis();

    double operator ()(int i, int j) const { auto it = rows[i].find(j); return it == rows[i].end()?0.0:it->second; }

    std::pair<int,int> findMarkowitzPivot();

    void pivot(int i,int j);

    int sparsestColInRow(int i);
    int sparsestRowInCol(int j);

    void addRowSparsityEntry(int i);
    void removeRowSparsityEntry(int i);
    void addColSparsityEntry(int j);
    void removeColSparsityEntry(int j);

    void setColBasic(int j, int i);
    void inactivateRow(int i);


    double meanColumnL0Norm() const;
    double meanColumnL1Norm() const;

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & cols & rows & basis & nAuxiliaryVars & F & L & U;
    }
};

template<class T>
std::ostream &operator <<(std::ostream &out, const TableauNormMinimiser<T> &tableau) {
    // print tableau
    for(int i=0; i<tableau.rows.size(); ++i) {
        out << tableau.L[i] << "\t<=\t";
        if(!tableau.isEqualityConstraint(i)) {
            out << "x(" << -tableau.basis[i] << ")";
        } else {
            out << "0   ";
        }
        out << "\t=\t";
        for(int j=0; j<tableau.cols.size(); ++j) {
            out << tableau(i,j) << "\t";
        }
        out << "+\t" << tableau.F[i] << " <= " << tableau.U[i] << std::endl;
    }
    // mark basic columns
    out << " \t \t \t   \t \t";
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
    L.reserve(problem.size());
    U.reserve(problem.size());
    int maxCol = 0; // maximum column id seen so far
    nAuxiliaryVars = 0;
    int constraintIndex = 0;
    for(const Constraint<T> &constraint: problem) {
        bool isActive = (constraint.upperBound == constraint.lowerBound);
        rows.emplace_back(constraint.coefficients, isActive);
        if (isActive) {
            addRowSparsityEntry(rows.size()-1);
            basis.push_back(INT_MAX); // fixed var
            F.push_back(-constraint.lowerBound); // transform all equality constraints to form DX - E = 0
            L.push_back(0);
            U.push_back(0);
        } else {
            basis.push_back(-(++nAuxiliaryVars));
            F.push_back(0);
            L.push_back(constraint.lowerBound);
            U.push_back(constraint.upperBound);
        }
        assert(rows.back().size() > 0);
        int rowMaxCol = rows.back().rbegin()->first;
        if(rowMaxCol > maxCol) maxCol = rowMaxCol;
        ++constraintIndex;
    }

    cols.resize(maxCol +1);
    // Now fill column information from rows
    for(int i=0; i<rows.size(); ++i) {
        for(auto [j, v]: rows[i]) {
            cols[j].insert(i);
        }
    }
    for(int j=0; j<=maxCol; ++j) {
        addColSparsityEntry(j);
    }
    std::cout << "Initial mean L0 norm = " << meanColumnL0Norm() << std::endl;
    std::cout << "Initial L1 norm = " << meanColumnL1Norm() << std::endl;
    findMinimalBasis();
}

template<class T>
void TableauNormMinimiser<T>::findMinimalBasis() {
    while(!rowsBySparsity.empty()) {
        auto [i,j] = findMarkowitzPivot();
        pivot(i,j);
    }
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
    T Mpipj = rows[pi].at(pj);

    // remove sparsity entries of all affected rows and columns
    for(auto [j, v] : rows[pi]) removeColSparsityEntry(j);

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

    // re-add previously removed sparsity entries
    for(int i : cols[pj]) {
        if(i != pi) {
            rows[i].erase(pj);
            if(rows[i].isActive) addRowSparsityEntry(i);
        }
    }
    setColBasic(pj, pi);
    for(const auto [j, Mpij]: rows[pi]) {
        if(j != pj) {
            cols[j].erase(pi);  // cols[] only stores non-reduced rows
            addColSparsityEntry(j);
        }
    }
    inactivateRow(pi);
}


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
        if(!cols[j].isBasic) {
            ++nNonBasic;
            for (int i: cols[j]) normSum += fabs(rows[i].at(j));
        }
    }
    return normSum / nNonBasic;
}


#endif //GLPKTEST_TABLEAUNORMMINIMISER_H
