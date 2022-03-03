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

// Takes an LP Problem and finds the basis that minimises the L1-norm of the tableau
// discounting any fixed-value variables
//  - all fixed-value variables non-basic?
//  - enforce integer coefficients?
// At a minimum, all valid pivots increase the matrix norm. So, find a pivot
// that decreases norm, pivot, repeat until no pivots exist.
// If we search pivots by row, and choose the pivot that minimises the norm
// then we can easily enforce the pivoting-out of fixed-value variables.
class TableauNormMinimiser {
public:
    static constexpr int maxVectorsToTest = 10;

    class Column: public std::set<int> {
    public:
        int sparseSize;
        std::list<int>::iterator sparsityEntry;
        bool isBasic;

        Column(const glp::SparseVec &col);
    };

    class Row: public std::map<int,double> {
    public:
        int sparseSize;
        std::list<int>::iterator sparsityEntry;
        bool isActive;

        Row(const glp::SparseVec &row, bool isActive);
    };

    std::vector<Column>     cols;
    std::vector<Row>        rows;

    std::vector<std::list<int>> colsBySparsity;
    std::vector<std::list<int>> rowsBySparsity;
    std::vector<int>            minimalBasis;

    TableauNormMinimiser(glp::Problem &problem);

     void findMinimalBasis();

    std::pair<int,int> findMarkowitzPivot();

    void pivot(int i,int j);

    int sparsestColInRow(int i);
    int sparsestRowInCol(int j);

    void updateRowSparsity(int i);
    void addRowSparsityEntry(int i);
    void removeRowSparsityEntry(int i);
    void updateColSparsity(int j);
    void addColSparsityEntry(int j);
    void removeColSparsityEntry(int j);

    void setColBasic(int j, int i);
    void inactivateRow(int i);

    double meanColumnNorm();
    double meanColumnL1Norm();
};


#endif //GLPKTEST_TABLEAUNORMMINIMISER_H
