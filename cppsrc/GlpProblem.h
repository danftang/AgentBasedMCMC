//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_GLPPROBLEM_H
#define GLPKTEST_GLPPROBLEM_H


#include <glpk.h>
#include "SparseVec.h"

class GlpProblem {
public:
    glp_prob *lp;
protected:
    SparseVec tmpStore;

public:
    GlpProblem(glp_prob *prob): tmpStore(glp_get_num_cols(prob)) {
        lp = prob;
    }


    int nConstraints()  { return glp_get_num_rows(lp); }
    int nVars()         { return glp_get_num_cols(lp); }

    const SparseVec &getObjective();
    const SparseVec &row(int i);
    const SparseVec &col(int j);
    double rowLowerBound(int i);
    double rowUpperBound(int i);
    double colLowerBound(int j);
    double colUpperBound(int j);

protected:
    SparseVec &getStore(int dim);
    SparseVec &getColStore()        { return getStore(nConstraints()); }
    SparseVec &getRowStore()        { return getStore(nVars()); }

};

std::ostream &operator <<(std::ostream &out, GlpProblem &prob);

#endif //GLPKTEST_GLPPROBLEM_H
