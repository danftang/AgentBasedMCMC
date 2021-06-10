//
// Created by daniel on 08/06/2021.
//
// Holds information regarding transition probabilities of
// a pivots on a column during the construction
// of an MCMC proposal.
//
//


#ifndef GLPKTEST_PROBABILISTICCOLUMNPIVOT_H
#define GLPKTEST_PROBABILISTICCOLUMNPIVOT_H


#include "ProposalPivot.h"

class ProbabilisticColumnPivot: public ProposalPivot {
public:
    static constexpr double tol = 1e-8;
    static constexpr double kappa = -0.001; // exponential coefficient for probabilities of choosing row based on change in feasibility

    bool leavingVarToUpperBound;
//    double transitionProb;
    std::vector<int> activeRows;
    std::vector<double> pivotPMF; // index is (2*activeRowIndex + toUpperBound), value is probability mass

    ProbabilisticColumnPivot(glp::Simplex &simplex, int j): ProbabilisticColumnPivot(simplex, j, simplex.tableauCol(j)) { }
    ProbabilisticColumnPivot(glp::Simplex &simplex, int j, std::vector<double> column); // j = column to pivot on

protected:

    static double iFeasibilityGradient(glp::Simplex &lp, int i, bool forward);
    double colFeasibilityGradient(glp::Simplex &lp, bool forward);
};


#endif //GLPKTEST_PROBABILISTICCOLUMNPIVOT_H
