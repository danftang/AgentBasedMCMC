//
// Created by daniel on 09/04/2021.
//

#include <cstdio>
#include <iostream>
#include <iomanip>
#include "GlpProblem.h"
#include "SparseVec.h"



const SparseVec &GlpProblem::getObjective() {
    SparseVec &obj = getRowStore();
    int nVars = glp_get_num_cols(lp);
    double c;
    for(int j=1; j<=nVars; ++j) {
        c = glp_get_obj_coef(lp, j);
        if(c != 0.0) obj.add(j-1, c);
    }
    return obj;
}


SparseVec &GlpProblem::getStore(int dim) {
    tmpStore.clear();
    tmpStore.n = dim;
    return tmpStore;
}

const SparseVec &GlpProblem::row(int i) {
    SparseVec &row = getRowStore();
    row.nnz = glp_get_mat_row(lp, i+1, row.ind, row.vec);
    return row;
}

const SparseVec &GlpProblem::col(int j) {
    SparseVec &col = getColStore();
    col.nnz = glp_get_mat_col(lp, j+1, col.ind, col.vec);
    return col;
}

double GlpProblem::rowLowerBound(int i) {
    return glp_get_row_lb(lp, i + 1);
}

double GlpProblem::rowUpperBound(int i) {
    return glp_get_row_ub(lp, i + 1);
}

double GlpProblem::colLowerBound(int j) {
    return glp_get_col_lb(lp, j + 1);
}

double GlpProblem::colUpperBound(int j) {
    return glp_get_col_ub(lp, j + 1);
}


std::ostream &operator <<(std::ostream &out, GlpProblem &prob) {
    switch(glp_get_obj_dir(prob.lp)) {
        case GLP_MIN: out << "Minimise "; break;
        case GLP_MAX: out << "Maximise "; break;
        default: out << "Unknown objective type ";
    }
    out << prob.getObjective() << std::endl;
    out << "Subject to:" << std::endl;
    for(int i=0; i<prob.nConstraints(); ++i) {
        out << prob.rowLowerBound(i) << " <= " << prob.row(i) << " <= " << prob.rowUpperBound(i) << std::endl;
    }
    for(int j=0; j<prob.nVars(); ++j) {
        out << prob.colLowerBound(j) << " <= x[" << j << "] <= " << prob.colUpperBound(j) << std::endl;
    }
    return out;
}

