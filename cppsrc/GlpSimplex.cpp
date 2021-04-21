//
// Created by daniel on 15/04/2021.
//

#include <iomanip>
#include "GlpSimplex.h"
extern "C" {
    #include "spxprob.h"
};


GlpSimplex::GlpSimplex(glp_prob *P) {
    spx_init_lp(this, P, 1);
    spx_alloc_lp(this);
    map = new int[1+P->m+P->n];
    spx_build_lp(this, P, 1, 1, map);
    spx_build_basis(this, P, map);
    pi = new double[m+1];
    pi[0] = 0.0;    // indicates not yet evaluated
    // b[1] = std::numeric_limits<double>::quiet_NaN();
}


GlpSimplex::~GlpSimplex() {
    spx_free_lp(this);
    delete [] map;
    delete [] pi;
}

void GlpSimplex::tableauRow(int i, double *trow) {
    double rho[m+1];
    spx_eval_rho(this, i, rho);
    spx_eval_trow(this, rho, trow);
}

double GlpSimplex::reducedObjective(int j) {
    if(pi[0] == 0.0) {
        spx_eval_pi(this, pi);
        pi[0] = 1.0;
    }
    return spx_eval_dj(this, pi, j);
}

// i,j is the pivot element, 1-based, (the j'th col of N replaces the i'th col of B)
// TODO: Better to use k coords on j?
// if leavingVarToUpperBound is true then the leaving var is set to its upper bound
// pivotCol - if non-null should be set to the current ftran of the incoming column,
// if null this will be calculated
void GlpSimplex::pivot(int i, int j, bool leavingVarToUpperBound, double *pivotCol) {
    const int boundFlag = leavingVarToUpperBound?1:0;
    double *pivCol;
    if(pivotCol == NULL) {
        pivCol = new double[m];
        tableauCol(j, pivCol);
    } else {
        pivCol = pivotCol;
    }
    spx_update_beta(this, b, i, boundFlag, j, pivCol);
    spx_update_invb(this,i,head[j+m]);
    spx_change_basis(this,i, boundFlag, j);
    pi[0]=0.0;
    if(pivotCol == NULL) delete [] pivCol;
}

void GlpSimplex::tableauCol(int j, double *col) {
    spx_eval_tcol(this, j, col);
}

std::ostream &operator <<(std::ostream &out, GlpSimplex &simplex) {
    const int nCols = simplex.n-simplex.m;
    double row[nCols+1];

    out << std::setprecision(5);

    // column labels
    out << std::setw(24) << "\t";
    for(int j=1; j<=nCols; ++j) {
        out << std::setw(12) << "x" << simplex.head[j+simplex.m]<<"\t";
    }
    out << std::setw(12) << "B" << std::endl;

    // col upper limits
    out << std::setw(24) << "\t";
    for(int j=1; j<=nCols; ++j) {
        out << std::setw(12) << simplex.u[simplex.head[j+simplex.m]] <<"\t";
    }
    out << std::endl;

    // cols on upper limits
    out << std::setw(24) << "\t";
    for(int j=1; j<=nCols; ++j) {
        out << std::setw(12) << (simplex.flag[j]==0?"\t":"--------")<<"\t";
    }
    out << std::endl;

    for(int i=1; i<=simplex.m; ++i) {
        out << std::setw(12) << simplex.l[simplex.head[i]] << " <= x" << simplex.head[i] << " = \t";
        simplex.tableauRow(i, row);
        for(int j=1; j<=nCols; ++j) {
            out << std::setw(12) << row[j] << "\t";
        }
        out << std::setw(12) << simplex.b[i] << " <= " << simplex.u[simplex.head[i]] << std::endl;
    }

    // reduced objective
    out << std::endl << std::setw(24) << "\t";
    for(int j=1; j<=nCols; ++j) {
        out << std::setw(12) << simplex.reducedObjective(j) <<"\t";
    }
    out << std::endl;

    // col lower limits
    out << std::endl << std::setw(24) << "\t";
    for(int j=1; j<=nCols; ++j) {
        out << std::setw(12) << simplex.l[simplex.head[j+simplex.m]] <<"\t";
    }
    out << std::endl;

    return out;
}
