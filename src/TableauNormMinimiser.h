//
// Created by daniel on 10/12/2021.
//

#ifndef GLPKTEST_TABLEAUNORMMINIMISER_H
#define GLPKTEST_TABLEAUNORMMINIMISER_H

#include <set>
#include <vector>
#include <map>
#include <list>
#include "glpkpp.h"
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
// M(X A)^T = F,  0 <= (X A) <= (H U-L)
//
// Given this, we can transform into the form
// Q X_N - F = (X A), 0 <= (X A) <= (H U-L)
// where A' is the non-fixed auxiliary variables, by taking
// columns of M over to the other side and adding singleton rows to M.
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

        explicit Column(const glp::SparseVec &col);
        Column(): isBasic(false) {};
    };

    class Row: public std::map<int,double> { // map from non-zero col index to element value
    public:
        std::list<int>::iterator sparsityEntry; // iterator into rowsBySaprsity entry
        bool isActive;

        Row(const glp::SparseVec &row, bool isActive);
    };

    std::vector<Column>     cols;   // tableau columns
    std::vector<Row>        rows;   // tableau rows

    std::vector<std::list<int>> colsBySparsity;     // only contains non-basic cols
    std::vector<std::list<int>> rowsBySparsity;     // only contains active rows
    std::vector<int>            minimalBasis;       // basic variables after minimisation: -ve is auxiliary by row index-1, otherwise col index

    TableauNormMinimiser(glp::Problem &problem);
    TableauNormMinimiser(ConvexPolyhedron &problem);

     void findMinimalBasis();

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

    double meanColumnL0Norm();
    double meanColumnL1Norm();
};


#endif //GLPKTEST_TABLEAUNORMMINIMISER_H
