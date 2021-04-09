//
// Created by daniel on 08/04/2021.
//
#include "math.h"
#include "glpkExtensions.h"

void glp_bfd_update(glp_prob *P, int j, int nnz, int *indices, double *vals) {
    bfd_update(P->bfd, j, nnz, indices, vals);
}
