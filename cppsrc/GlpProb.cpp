//
// Created by daniel on 09/04/2021.
//

#include <cstdio>
#include <iostream>
#include <iomanip>
#include "GlpProb.h"
#include "SparseVec.h"

SparseVec &GlpProb::getTableauRow(int i) {
    int k = glp_get_bhead(lp, i);
    rowVec.sparseSize = glp_eval_tab_row(lp, k, rowVec.indices, rowVec.vals);
    rowVec.add(k, 1.0);
    return rowVec;
}


void GlpProb::printTableau() {
    int nRows = glp_get_num_rows(lp);
    int nCols = glp_get_num_cols(lp) + nRows;
    int i, j, k;
    double Mij[nCols - nRows];

    std::cout << std::setprecision(5);

    std::cout << "Tableau: " << nRows << " x " << nCols << std::endl;
    for (i = 1; i <= nRows; i++) {
        std::cout << std::setw(12) << glp_get_row_ub(lp, i) << "\t";
    }
    for (j = 1; j <= nCols - nRows; j++) {
        std::cout << std::setw(12) << glp_get_col_ub(lp, j) << "\t";
    }
    std::cout << std::endl;
    for (i = 1; i <= nRows; i++) {
        switch (glp_get_row_stat(lp, i)) {
            case GLP_BS:
                std::cout <<"************\t";
                break;
            case GLP_NU:
                std::cout <<"------------\t";
                break;
            default:
                std::cout <<"             \t";
        }
    }
    for (j = 1; j <= nCols - nRows; j++) {
        switch (glp_get_col_stat(lp, j)) {
            case GLP_BS:
                std::cout <<"************\t";
                break;
            case GLP_NU:
                std::cout <<"------------\t";
                break;
            default:
                std::cout <<"             \t";
        }
    }
    std::cout << std::endl;
    for (i = 1; i <= nRows; i++) {
        std::cout << getTableauRow(i) ;
        k = glp_get_bhead(lp, i);
        if (k <= nRows) {
            std::cout << glp_get_row_prim(lp, k) << std::endl;
        } else {
            std::cout << glp_get_col_prim(lp, k - nRows) << std::endl;
        }
    }
    for (i = 1; i <= nRows; i++) {
        std::cout << std::setw(12) << glp_get_row_dual(lp, i) << "\t";
    }
    for (j = 1; j <= nCols - nRows; j++) {
        std::cout << std::setw(12) << glp_get_col_dual(lp, j) << "\t";
    }
    std::cout << glp_get_obj_val(lp) << std::endl;
    for (i = 1; i <= nRows; i++) {
        if (glp_get_row_stat(lp, i) == GLP_NL) {
            std::cout <<"------------\t";
        } else {
            std::cout <<"            \t";
        }
    }
    for (j = 1; j <= nCols - nRows; j++) {
        if (glp_get_col_stat(lp, j) == GLP_NL) {
            std::cout <<"------------\t";
        } else {
            std::cout <<"            \t";
        }
    }
    std::cout << std::endl;
    for (i = 1; i <= nRows; i++) {
        std::cout << std::setw(12) << glp_get_row_lb(lp, i) << "\t";
    }
    for (j = 1; j <= nCols - nRows; j++) {
        std::cout << std::setw(12) << glp_get_col_lb(lp, j) << "\t";
    }
    std::cout << std::endl;
}