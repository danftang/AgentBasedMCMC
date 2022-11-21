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
#include "EqualityConstraints.h"
#include "ConstrainedFactorisedDistribution.h"
#include "include/StlStream.h"

// Takes an LP Problem and finds the basis that pivots out as many
// fixed-variables as possible, while attempting to keep the L0-norm of the tableau
// as low as possible.
//
// Starting with a set of linear constraints on X
//  DX = F
//
// If D has m rows, this class chooses m linearly independent columns
// of D and pivots on rows in D so that X can be split into X_N and
// X_B so that
//
// D'X = (I  M)(X_B) = F'
//             (X_N)
//
// All fixed variables are assumed to be functions of integers only.
// Pivoting only occurs on elements that have value +-1 to retain
// integer coefficients.
//
// We want to pivot in
// such a way as to maintain the sparsity of these matrices. We do this using a
// greedy algorithm based on the Markovitz criterion. See
//  Maros, I., 2002, "Computational techniques of the sipmlex method", Springer
//
//
template<class COEFF>
class TableauNormMinimiser {
public:
    template<class DOMAIN>
    static void factorise(ConstrainedFactorisedDistribution<DOMAIN,COEFF> &distribution) {
        TableauNormMinimiser minimiser(distribution);
        distribution.constraints = minimiser.getFactorisedConstraints();
    }

public:

    static constexpr int MAX_VECTORS_TO_TEST = 15;
    static constexpr int UNREDUCED = -1;


    class Column: public std::set<int> { // set of row indices of non-zero entries, +ve is equality constraint, -ve is factor dependency
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
        }
    };

    class Row: public std::map<int,COEFF> { // map from non-zero col index to element value
    public:
        std::list<int>::iterator sparsityEntry; // iterator into rowsBySaprsity entry
        bool isActive;                          // true if this row is an unreduced equality constraint

        Row(const SparseVec<COEFF> &row, bool isActive);
        Row() {}

        Row &operator *=(COEFF c);

        // dot product
        template<typename OTHER, typename E = decltype(std::declval<COEFF &>() = std::declval<COEFF>() * std::declval<OTHER>()[0])>
        COEFF operator *(const OTHER &other) const {
            COEFF dotProd = 0;
            for(const std::pair<int,COEFF> &element: *this) {
                dotProd += element.second * other[element.first];
            }
            return dotProd;
        }

    private:
        friend class boost::serialization::access;

        template <typename Archive>
        void serialize(Archive &ar, const unsigned int version) {
            ar & static_cast<std::map<int,COEFF> &>(*this) & isActive;
        }
    };

    std::vector<Column>     cols;   // tableau columns
    std::vector<Row>        rows;   // tableau rows

    std::vector<std::list<int>> colsBySparsity;     // cols sorted by sparsity, only contains non-basic cols
    std::vector<std::list<int>> rowsBySparsity;     // rows sorted by sparsity, only contains active rows
    std::vector<int>            basicVars;            // col index of basic var by row. -1 means unreduced
    std::vector<COEFF>              F;                  // constant in linear equations by row (see intro above)

    template<class DOMAIN>
    TableauNormMinimiser(const ConstrainedFactorisedDistribution<DOMAIN,COEFF> &distribution);

    TableauNormMinimiser()=default;

    void findMinimalBasis();
    std::vector<SparseVec<COEFF>> getBasisVectors();
    EqualityConstraints<COEFF> getFactorisedConstraints();
//    template<class DOMAIN>
//    Basis<DOMAIN> getBasis();
    SparseVec<COEFF> getOrigin();

    double operator ()(int i, int j) const { auto it = rows[i].find(j); return it == rows[i].end()?0.0:it->second; }

    template<typename VEC>
    void snapToSubspace(VEC &X) const;

protected:


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


};

template<class T>
std::ostream &operator <<(std::ostream &out, const TableauNormMinimiser<T> &tableau) {
    // print tableau
    for(int i=0; i<tableau.rows.size(); ++i) {
        out << "x(" << tableau.basicVars[i] << ")";
        out << "\t=\t";
        for(int j=0; j<tableau.cols.size(); ++j) {
            out << tableau(i,j) << "\t";
        }
        out << "+\t" << tableau.F[i] << std::endl;
    }
    // mark basic columns
    out << "\t   \t \t";
    for(int j=0; j<tableau.cols.size(); ++j) {
        out << (tableau.cols[j].isBasic?"B\t":"-\t");
    }
    out << std::endl;
    out << tableau.cols << std::endl;
    return out;
}

// Turns the supplied distribution into a set of equality constraints
// and factor dependencies.
// -ve entries in columns are factor dependencies and non-negative entries
// are constraints.
template<class T>
template<class D>
TableauNormMinimiser<T>::TableauNormMinimiser(const ConstrainedFactorisedDistribution<D,T> &distribution)
{
    basicVars.reserve(distribution.constraints.size());
    F.reserve(distribution.constraints.size());
    rows.reserve(distribution.constraints.size());
    // initialise row information
    for(const EqualityConstraint<T> &constraint: distribution.constraints) {
        rows.emplace_back(constraint.coefficients, true);
        addRowSparsityEntry(rows.size()-1);
        basicVars.push_back(UNREDUCED); // unreduced row
        F.push_back(-constraint.constant); // transform all equality constraints to form DX + F = 0
        assert(rows.back().size() > 0);
    }

    // Now fill column information from rows and factor dependencies
    cols.resize(D::dimension);
    for (int factorIndex = 0; factorIndex < distribution.factors.size(); ++factorIndex) {
            for (int basisIndex: distribution.factors[factorIndex].dependencies) {
                cols[basisIndex].insert(-1 - factorIndex);
            }
    }

    for(int i=0; i<rows.size(); ++i) {
        for(auto [j, v]: rows[i]) {
            cols[j].insert(i);
        }
    }
//    for(int j=0; j<cols.size(); ++j) {
//        addColSparsityEntry(j);
//    }
//    for(int j=cols.size()-1; j>=0; --j) {
//        addColSparsityEntry(j);
//    }

    std::vector<int> colIndices(cols.size());
    std::generate(colIndices.begin(), colIndices.end(), [i = 0]() mutable { return i++; });
    std::shuffle(colIndices.begin(), colIndices.end(), Random::gen);
    for(int j: colIndices) {
        addColSparsityEntry(j);
    }



    std::cout << "Initial mean column L0 norm = " << meanColumnL0Norm() << std::endl;
    std::cout << "Initial column L1 norm = " << meanColumnL1Norm() << std::endl;
    findMinimalBasis();
}

template<class T>
void TableauNormMinimiser<T>::findMinimalBasis() {
    std::cout << "Finding minimal basis for tableau of size " << rows.size() << " x " << cols.size() << std::endl;
//    std::cout << colsBySparsity << std::endl << std::endl;
//    std::cout << rowsBySparsity << std::endl << std::endl;

    while(!rowsBySparsity.empty()) {
        auto [i,j] = findMarkowitzPivot();
        pivot(i,j);
    }

    std::vector<int> basics = basicVars;
    std::sort(basics.begin(), basics.end());
//    std::cout << "Basic vars: " << basics << std::endl;
    std::cout << "Reduced mean column L0 norm = " << meanColumnL0Norm() << std::endl;
    std::cout << "Reduced column L1 norm = " << meanColumnL1Norm() << std::endl;

}

template<class T>
std::pair<int,int> TableauNormMinimiser<T>::findMarkowitzPivot() {
    int k=1;
    int lowestScore = INT32_MAX;
    int vectorsTested = 0;
    std::pair<int,int> bestPivot(-1,-1);
    while(k<colsBySparsity.size() && k<rowsBySparsity.size()) {
//        std::vector<int> sparseCols(colsBySparsity[k].begin(), colsBySparsity[k].end());
//        std::shuffle(sparseCols.begin(), sparseCols.end(), Random::gen);
//        std::vector<int> sparseRows(rowsBySparsity[k].begin(), rowsBySparsity[k].end());
//        std::shuffle(sparseRows.begin(), sparseRows.end(), Random::gen);
        for(int i: rowsBySparsity[k]) {
            int bestCol = sparsestColInRow(i); // TODO: if k=1 or 2 best col is the densest
            if(bestCol != -1) {
                int score = (rows[i].size() - 2) * (cols[bestCol].size() - 1);
                if (score < lowestScore) {
                    bestPivot = std::pair<int, int>(i, bestCol);
                    if (score <= (k - 2) * (k - 1)) {
//                        std::cout << "Found best pivot on row of density " << k << std::endl;
                        return bestPivot;
                    }
                    lowestScore = score;
                }
                if (++vectorsTested >= MAX_VECTORS_TO_TEST) {
//                    std::cout << "Found pivot at " << bestPivot << " of density " << rows[bestPivot.first].size() << " x " << cols[bestPivot.second].size() << std::endl;
                    return bestPivot;
                }
            }
        }

        for(int j: colsBySparsity[k]) {
            int bestRow = sparsestRowInCol(j);
            if(bestRow != -1) {
                int score = (rows[bestRow].size() - 2) * (cols[j].size() - 1);
                if (score < lowestScore) {
                    bestPivot = std::pair<int, int>(bestRow, j);
                    if (score <= (k - 1) * (k - 1)) {
//                        std::cout << "Found best pivot on col " << j << " of density " << k << std::endl;
                        return bestPivot;
                    }
                    lowestScore = score;
                }
                if (++vectorsTested >= MAX_VECTORS_TO_TEST) {
//                    std::cout << "Found pivot at " << bestPivot << " of density " << rows[bestPivot.first].size() << " x " << cols[bestPivot.second].size() << std::endl;
                    return bestPivot;
                }
            }
        }
        ++k;
    }
//    std::cout << "Found pivot at " << bestPivot << " of density " << rows[bestPivot.first].size() << " x " << cols[bestPivot.second].size() << std::endl;
    return bestPivot;
}


template<class T>
void TableauNormMinimiser<T>::pivot(int pi, int pj) {
    assert(rows[pi].isActive);
    T Mpipj = rows[pi].at(pj);

    // remove sparsity entries of all affected columns
    for(const auto &[j, v] : rows[pi]) removeColSparsityEntry(j);

    rows[pi] *= -1.0/Mpipj;     // divide row through to make the pivot point -1
    F[pi] *= -1.0/Mpipj;

    for (int i: cols[pj]) {
        if (i != pi) {
            if (i >= 0) { // row is equality constraint, so do linear pivot
                if(rows[i].isActive) removeRowSparsityEntry(i);
                T Mipj = rows[i].at(pj);
                for (const auto &[j, Mpij]: rows[pi]) {
                    if (j != pj) {
                        auto MijIt = rows[i].find(j);
                        if (MijIt == rows[i].end()) { // fill-in
                            rows[i][j] = Mipj * Mpij;
                            cols[j].insert(i);
                        } else {
                            T newMij = (MijIt->second += Mipj * Mpij);
                            if (newMij == 0) { // drop-out
                                cols[j].erase(i);
                                rows[i].erase(j);
                            }
                        }
                    } else {
                        // drop out without testing as we're on the pivot col
                        // don't erase col entry yet as we're iterating through it
                        rows[i].erase(j);
                    }
                }
                if (rows[i].isActive) addRowSparsityEntry(i);
                F[i] += Mipj * F[pi];
            } else { // row is factor dependency so do union pivot
                for (const auto &[j, Mpij]: rows[pi]) cols[j].insert(i);
            }
        }
    }

//    // re-add previously removed sparsity entries
//    for(int i : cols[pj]) {
//        if(i != pi) {
//            rows[i].erase(pj);
//            if(rows[i].isActive) addRowSparsityEntry(i);
//        }
//    }
    setColBasic(pj, pi); // removes entries from col and records in basis
    // add new column sparsity entries
    for(const auto [j, Mpij]: rows[pi]) {
        if(j != pj) {
//            cols[j].erase(pi);  // cols[] only stores non-reduced rows (triangular factorisation)
            addColSparsityEntry(j);
        }
    }
    inactivateRow(pi); // sets inactive and removes row sparsity entry
}


template<class T>
void TableauNormMinimiser<T>::addRowSparsityEntry(int i) {
    assert(i>=0);
    Row &row = rows[i];
    if(rowsBySparsity.size() <= row.size()) rowsBySparsity.resize(row.size()+1);
    rowsBySparsity[row.size()].push_front(i);
    row.sparsityEntry = rowsBySparsity[row.size()].begin();
}

template<class T>
void TableauNormMinimiser<T>::removeRowSparsityEntry(int i) {
    assert(i>=0);
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
    basicVars[i] = j;
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
        if(colSize < minColSparsity && !cols[j].isBasic && Mij*Mij == 1) {
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
        if(i >=0) {
            int rowSize = rows[i].size();
            if (rowSize < minRowSparsity && rows[i].isActive && fabs((double) rows[i].at(j)) == 1.0) {
                minRowSparsity = rowSize;
                mini = i;
            }
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
            for (int i: cols[j]) {
                if(i>=0) normSum += fabs((double)rows[i].at(j));
            }
        }
    }
    return normSum / nNonBasic;
}

template<class T>
std::vector<SparseVec<T>> TableauNormMinimiser<T>::getBasisVectors() {
//    std::cout << "Starting with tableaux: " << std::endl << *this << std::endl;
//    std::cout << rowsBySparsity << std::endl;
//    std::cout << colsBySparsity << std::endl;
//    std::cout << cols << std::endl;
//    findMinimalBasis();
//    std::cout << "Finished with tableaux: " << std::endl << *this << std::endl;
//    std::cout << rowsBySparsity << std::endl;
//    std::cout << colsBySparsity << std::endl;
//    std::cout << cols << std::endl;

    std::vector<SparseVec<T>> basisVectors;
    basisVectors.reserve(cols.size());
    for(int j=0; j<cols.size(); ++j) {
        if(!cols[j].isBasic) {
            SparseVec<T> &newBasis = basisVectors.emplace_back();
            newBasis.insert(j, 1); // element from the identity (insert first to ensure canonical)
            for (int i: cols[j]) {
                if (i >= 0) {
                    assert(basicVars[i] != UNREDUCED); // every row should have a basic variable
                    newBasis.insert(basicVars[i], rows[i][j]);
                }
            }
        }
    }
    return basisVectors;
}

//template<class T>
//template<class DOMAIN>
//Basis<DOMAIN> TableauNormMinimiser<T>::getBasis() {
//    std::vector<SparseVec<T>> basisVecs = getBasisVectors(DOMAIN::dimension);
//    DOMAIN origin; // assumes default constructor generates zero vector
//    for(int row =0; row < basicVars.size(); ++row) origin[basicVars[row]] = F[row];
//    return Basis(basisVecs, origin);
//}


// returns the value of the domain when all non-basic variables are zero
template<class T>
SparseVec<T> TableauNormMinimiser<T>::getOrigin() {
    SparseVec<T> origin;
    for(int row =0; row < basicVars.size(); ++row) origin.insert(basicVars[row], F[row]);
    return origin;
}


// Get the solution as a set of constraints, with each constraint having a basic variable
// that is not present in any other constraint
template<class T>
EqualityConstraints<T> TableauNormMinimiser<T>::getFactorisedConstraints() {
    EqualityConstraints<T> constraints;
    constraints.reserve(rows.size());
    for(int i = 0; i<rows.size(); ++i) constraints.emplace_back(rows[i], F[i]);
    return constraints;
}


// Updates the basic variables of X, leaving the non-basic variables unchanged
// so that X satisfies all constraints.
// TODO: should this be elsewhere?
template<class T>
template<typename VEC>
void TableauNormMinimiser<T>::snapToSubspace(VEC &X) const {
    for(int i=0; i<rows.size(); ++i) {
        int basisIndex = basicVars[i];
        X[basisIndex] = 0;
        X[basisIndex] = rows[i]*X + F[i];
    }
}


#endif //GLPKTEST_TABLEAUNORMMINIMISER_H
