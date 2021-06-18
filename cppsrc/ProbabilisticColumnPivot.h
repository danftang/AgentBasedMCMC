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
    static constexpr double kappa = -8.0; // exponential coefficient for probabilities of choosing row based on change in feasibility

//    double transitionProb;
    glp::Simplex &simplex;
    std::vector<int> nonZeroRows;
    std::vector<double> pivotPMF; // index is (2*activeRowIndex + toUpperBound), value is probability mass

    ProbabilisticColumnPivot(glp::Simplex &simplex, int j): ProbabilisticColumnPivot(simplex, j, simplex.tableauCol(j)) { }
    ProbabilisticColumnPivot(glp::Simplex &simplex, int j, std::vector<double> column): ProposalPivot(-1, j, std::move(column)), simplex(simplex) {
        chooseRow();
    }
    ProbabilisticColumnPivot(glp::Simplex &simplex): simplex(simplex) {
        chooseCol();
        col = simplex.tableauCol(j);
        chooseRow();
    }

    // Test stuff
    double feasibility(double deltaj);
protected:

//    double iFeasibilityGradient(int i, bool forward);
//    double colFeasibilityGradient(bool forward);
    double colFeasibilityGradient(double deltaj);
    bool isActive(int pmfIndex);
    void chooseCol();
    void chooseRow();

    double feasibilityGradient(double v, double lowerBound, double upperBound);
};


#endif //GLPKTEST_PROBABILISTICCOLUMNPIVOT_H
