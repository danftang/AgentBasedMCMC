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
#include "SimplexMCMC.h"
#include "Phase1Pivot.h"

class ProbabilisticColumnPivot: public Phase1Pivot {
public:
    static constexpr double kappa = -8.0; // exponential coefficient for probabilities of choosing row based on change in infeasibility
    static constexpr double alpha = 100.0; // relative probability of choosing column with non-zero reduced objective


    explicit ProbabilisticColumnPivot(glp::Simplex &simplex);
//    ProbabilisticColumnPivot(glp::Simplex &simplex, int j): ProposalPivot(-1, j, simplex.tableauCol(j)), simplex(simplex) {
//        chooseRow();
//    }
//    ProbabilisticColumnPivot(glp::Simplex &simplex, int i, int j): ProposalPivot(i, j, simplex.tableauCol(j)), simplex(simplex) { }

    void chooseCol();
    void chooseRow();

    // Test stuff
    double infeasibility(double deltaj);
//    int infeasibilityCount(double deltaj);


protected:

    bool isActive(int pmfIndex);
    void calcAcceptanceContrib();
};


#endif //GLPKTEST_PROBABILISTICCOLUMNPIVOT_H
