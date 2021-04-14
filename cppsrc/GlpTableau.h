//
// Created by daniel on 12/04/2021.
//

#ifndef GLPKTEST_GLPTABLEAU_H
#define GLPKTEST_GLPTABLEAU_H

#include "glpk.h"
#include "SparseVec.h"

class GlpTableau {
public:
    glp_prob *lp;

public:
    GlpTableau(glp_prob *prob) { lp = prob; }

    int nRows()         { return glp_get_num_rows(lp); }
    int nCols()         { return glp_get_num_cols(lp) + glp_get_num_rows(lp); }

    void row(int i, SparseVec &result);    // row of the tableau
    void col(int k, SparseVec &result);    // col of the tableau
    void MCol(int k, SparseVec &result);         // col of the original problem in tmpStore (first nRows vars are artificial)
    void MCol(int k, double *result);

    void reducedObjective(SparseVec &result);
    bool isBasic(int j);
    int basicRow(int k);    // if k is basic col, returns non-zero row entry else 0
    int basicCol(int i)     { return glp_get_bhead(lp, i+1)-1; } // returns the basic col of row i
    int colState(int k);
    void setColState(int k, int state);
    double upperBound(int k);
    double lowerBound(int k);
    double objectiveVal()   { return glp_get_obj_val(lp); }
    double basicVal(int i);
    bool pivot(int i, int j);

    int colType(int k);

    void setBasicRow(int i, int k);

    double colVal(int k);

    void incrementBasicVal(int i, double increment);
};

std::ostream &operator <<(std::ostream &out, GlpTableau &tableau);

#endif //GLPKTEST_GLPTABLEAU_H
