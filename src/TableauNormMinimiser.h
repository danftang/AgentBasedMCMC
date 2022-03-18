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

    class Row: public std::map<int,double> { // map from non-zero col index to element value
    public:
        std::list<int>::iterator sparsityEntry; // iterator into rowsBySaprsity entry
        bool isActive;

        Row(const SparseVec<double> &row, bool isActive);

        Row &operator *=(double c);
    };

    std::vector<Column>     cols;   // tableau columns
    std::vector<Row>        rows;   // tableau rows

    std::vector<std::list<int>> colsBySparsity;     // only contains non-basic cols
    std::vector<std::list<int>> rowsBySparsity;     // only contains active rows
    std::vector<int>            basis;            // basic variables: -ve means auxiliary, otherwise col index
    int                         nAuxiliaryVars;     // number of auxiliary variables
    std::vector<double>         F;                  // aggregated value of all non-basic fixed vars
    std::vector<double>         Hc;                  // upper bounds of columns
    std::vector<double>         Ha;                  // upper bounds of auxiliaries

//    TableauNormMinimiser(glp::Problem &problem);
    TableauNormMinimiser(const ConvexPolyhedron &problem);

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

    friend std::ostream &operator <<(std::ostream &out, const TableauNormMinimiser &tableau);
};


#endif //GLPKTEST_TABLEAUNORMMINIMISER_H
