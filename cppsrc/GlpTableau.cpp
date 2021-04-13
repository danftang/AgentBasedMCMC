//
// Created by daniel on 12/04/2021.
//

#include <iostream>
#include <iomanip>
#include "GlpTableau.h"

extern "C" {
    #include "prob.h"
}

void GlpTableau::row(int i, SparseVec &result) {
    SparseVec c(nRows());
    result.clear();
    result.n = nCols();
    double Mik;

    // using btran
    // row = e_i B^-1 N = rho N
    // where rho = e_i B^-1
    double rho[nRows()];
    for(int j=0; j<nRows(); ++j) rho[j] = 0.0;
    rho[i] = 1.0;
    bfd_btran(lp->bfd, rho-1);
    for(int k=0; k<nCols(); ++k) {
        MCol(k, c);
        Mik = c.dotProd(rho);
        if(fabs(Mik) > 1e-15) result.add(k, Mik);
    }

    // using eval_tab_row
//    int k = glp_get_bhead(lp, i+1)-1;
//    result.nnz = glp_eval_tab_row(lp, k+1, result.ind, result.vec);
//    result.add(k, 1.0);
}

void GlpTableau::col(int k, SparseVec &result) {
    result.clear();
    result.n = nRows();

    // using bfd_ftran
    double dense[nRows()];
    MCol(k, result);
    bfd_ftran(lp->bfd, dense-1);
//    glp_ftran(lp, dense);
    for(int i=0; i<nRows(); ++i) {
        if(fabs(dense[i]) > 1e-15) result.add(i, dense[i]);
    }

    // using eval_tab_col
//    if(isBasic(k)) {
//        col.add(basicRow(k),1.0);
//    } else {
//        col.nnz = glp_eval_tab_col(lp, k + 1, col.ind, col.vec);
//        for(int i=1; i<=col.nnz; ++i) {
//            col.ind[i] = basicRow(col.ind[i]-1)+1;
//        }
//    }
}


void GlpTableau::reducedObjective(SparseVec &result) {
    for (int i = 0; i < nRows(); i++) {
        result.add(i, glp_get_row_dual(lp, i+1));
    }
    for (int j = nRows(); j < nCols(); j++) {
        result.add(j, glp_get_col_dual(lp, j-nRows()+1));
    }
}


bool GlpTableau::isBasic(int k) {
    return colState(k) == GLP_BS;
}


int GlpTableau::basicRow(int k) {
    if(k < nRows()) {
        return glp_get_row_bind(lp, k + 1)-1;
    }
    return glp_get_col_bind(lp, k - nRows() + 1)-1;
}


int GlpTableau::colState(int k) {
    if(k < nRows()) return glp_get_row_stat(lp,k+1);
    return glp_get_col_stat(lp,k-nRows()+1);
}

double GlpTableau::upperBound(int k) {
    if(k < nRows()) return glp_get_row_ub(lp, k+1);
    return glp_get_col_ub(lp, k-nRows()+1);
}

double GlpTableau::lowerBound(int k) {
    if(k < nRows()) return glp_get_row_lb(lp, k+1);
    return glp_get_col_lb(lp, k-nRows()+1);
}

// i is row
double GlpTableau::basicVal(int i) {
    int k = glp_get_bhead(lp, i+1);
    if (k <= nRows()) return glp_get_row_prim(lp, k);
    return glp_get_col_prim(lp, k - nRows());
}

// Removes the basic variable associated with row i from the
// basis and adds k. The entering variable is moved away from its
// bound until the leaving variable hits one of its bounds.
// However, if the leaving variable has only one bound, it
// is set to that bound.
// The entering variable, k, becomes basic.
bool GlpTableau::pivot(int i, int k) {
    int incomingVarState = colState(k);
    int leavingVar = basicCol(i);
    switch(incomingVarState) {
        case GLP_BS: if(leavingVar == k) return true; else return false;
        case GLP_NF: return false;
        case GLP_NS: return false;
    }
    col(k,tmpColStore());
    bool pivotPointIsPositive = tmpStore[i] > 0.0;
//    std::cout << "Pivoting row " << i << " to " << c << std::endl;
    int failed = bfd_update(lp->bfd, i+1, tmpStore.nnz, tmpStore.ind, tmpStore.vec);
    if(failed != 0) return false;

    // make leaving var non-basic
    if(incomingVarState == GLP_NL) {
        setColState(leavingVar, pivotPointIsPositive?GLP_NL:GLP_NU);
    } else {
        setColState(leavingVar, pivotPointIsPositive?GLP_NU:GLP_NL);
    }
    if(leavingVar < nRows()) {
        lp->row[leavingVar+1]->bind = 0;
    } else {
        lp->col[leavingVar-nRows()+1]->bind = 0;
    }

    // make entering variable basic
    setColState(k, GLP_BS);
    lp->head[i+1] = k+1;
    if(k < nRows()) {
        lp->row[k+1]->bind = i+1;
    } else {
        int j = k - nRows() + 1;
        lp->col[j]->bind = i+1;
    }

    lp->valid = 1;
    return true;
}


void GlpTableau::setColState(int k, int state) {
    if(k < nRows()) glp_set_row_stat(lp, k+1, state); else glp_set_col_stat(lp, k-nRows()+1, state);
}


void GlpTableau::MCol(int k, SparseVec &result) {
    result.clear();
    result.n = nRows();
    if(k < nRows()) {
        result.add(k,1.0);
    } else {
        result.nnz = glp_get_mat_col(lp, k - nRows() + 1, result.ind, result.vec);
    }
}

void GlpTableau::MCol(int k, double *result) {
    MCol(k, tmpColStore());
    tmpStore.toDense(result);
}

std::ostream &operator <<(std::ostream &out, GlpTableau &tableau) {
    SparseVec tmpSparse(tableau.nCols());
    int nRows = tableau.nRows();
    int nCols = tableau.nCols();
    int i, k;

    out << std::setprecision(5);

    out << "Tableau: " << nRows << " x " << nCols << std::endl;
    for (k = 0; k < nCols; k++) {
        out << std::setw(12) << tableau.upperBound(k) << "\t";
    }
    out << std::endl;
    for (k = 0; k < nCols; k++) {
        switch (tableau.colState(k)) {
            case GLP_BS:
                out <<"************\t";
                break;
            case GLP_NU:
                out <<"------------\t";
                break;
            default:
                out <<"             \t";
        }
    }
    out << std::endl;
    for (i = 0; i < nRows; i++) {
        tableau.row(i,tmpSparse);
        out << tmpSparse << tableau.basicVal(i) << std::endl;
    }
    tableau.reducedObjective(tmpSparse);
    out << std::endl << tmpSparse << tableau.objectiveVal() << std::endl;
    for (k = 0; k < nCols; k++) {
        if (tableau.colState(k) == GLP_NL) out << "------------\t"; else out << "            \t";
    }
    std::cout << std::endl;
    for (k = 0; k < nCols; k++) {
        out << std::setw(12) << tableau.lowerBound(k) << "\t";
    }
    out << std::endl;
    return out;
}
