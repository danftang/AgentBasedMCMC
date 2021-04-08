//
// Created by daniel on 08/04/2021.
//
#include "math.h"
#include "glpkExtensions.h"

void toDense(int nnz, const int idx[], const double vals[], int n, double dense[]) {
    for(int i=0; i<n; ++i) {
        dense[i] = 0.0;
    }
    for(int i=1; i<=nnz; ++i) {
        dense[idx[i]-1] = vals[i];
    }
}

void glp_print_tableau(glp_prob *lp) {
    int nnz;
    int nCols = glp_get_num_cols(lp);
    int nRows = glp_get_num_rows(lp);
    int idx[nCols];
    int k;
    double Mij[nCols];
    double row[nCols];

    printf("Tableau: %d %d \n", nRows, nCols);
    for(int i=1; i <= nRows; i++) {
        k = glp_get_bhead(lp, i);
        nnz = glp_eval_tab_row(lp, k, idx, Mij);
        toDense(nnz, idx, Mij, nCols, row);
        printf("%d : ", k);
        for(int j=0; j<nCols; j++) {
            printf("%f ",row[j]);
        }
        printf(" : %f\n", (k-nRows)>0?glp_get_col_prim(lp,  k-nRows):NAN);
    }
}

