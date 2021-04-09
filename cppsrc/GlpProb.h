//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_GLPPROB_H
#define GLPKTEST_GLPPROB_H


#include <glpk.h>
#include "SparseVec.h"

class GlpProb {
public:
    glp_prob *lp;
    SparseVec rowVec;
    SparseVec colVec;

    GlpProb(glp_prob *prob): rowVec(glp_get_num_rows(prob) + glp_get_num_cols(prob)), colVec(glp_get_num_rows(prob)) {
        lp = prob;
    }

    SparseVec &getTableauRow(int i);

    void printTableau();
};


#endif //GLPKTEST_GLPPROB_H
