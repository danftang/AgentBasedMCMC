//
// Created by daniel on 08/04/2021.
//

#ifndef GLPKTEST_GLPKEXTENSIONS_H
#define GLPKTEST_GLPKEXTENSIONS_H

#include "glpk.h"
#include "spxprob.h"

void glp_bfd_update(glp_prob *P, int j, int nnz, int *indices, double *vals);

#endif //GLPKTEST_GLPKEXTENSIONS_H



