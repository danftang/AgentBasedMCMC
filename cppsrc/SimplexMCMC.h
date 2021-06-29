//
// Created by daniel on 15/04/2021.
//

#ifndef GLPKTEST_SIMPLEXMCMC_H
#define GLPKTEST_SIMPLEXMCMC_H

#include <vector>
#include <set>
#include <map>
#include <random>
#include "glpkpp.h"
#include "ProposalPivot.h"

class SimplexMCMC: public glp::Simplex {
public:
    static constexpr double fractionalK = 0.1;
//    static constexpr double tol = 1e-8;

//    using glp::Simplex::pivot;

    std::vector<int>    colLastNonZero;             // row-index of the last non-zero entry in this col
    std::vector<int>    rowLatestCompletionPivot; // k-col-index of the earliest col of the final basis that can pivot on this row
    std::vector<int>    latestCompletionBegin;   // k-col-index of the earliest col of the final completion
    std::vector<double> lnRowPivotCount;            // ln of number of possible choices of next in-sequence var of degenerate state
   // std::set<int>       finalBasis;                 // the final degenerate state (ordered)
//    std::map<int,int>   orderedBasis;               // the current basis ordered map from k-index to i-index
    int nSamples = 0;
    int nRejections = 0;
    int fractionalRunLength = 0;
    std::function<double (const std::vector<double> &)> logProbFunc;

//    BasisProbability probability;


    SimplexMCMC(glp::Problem &prob, const std::function<double (const std::vector<double> &)> &logProb);

    double lnDegeneracyProb();
    double lnProb() { return logProbFunc(X()); }
    double lnFractionalPenalty();

    void nextSample();
//    double reverseTransitionProb(ProposalPivot proposal);
    void pivot(const ProposalPivot &piv) { this->glp::Simplex::pivot(piv.i, piv.j, piv.col, piv.leavingVarToUpperBound); }


//    std::vector<int> calcPivotRows(int j, const std::vector<double> &colVec);
//    std::vector<int> calcPivotRows(int j) { return calcPivotRows(j,tableauCol(j)); }

    void randomWalk();

    static glp::Problem &initialiseProblem(glp::Problem &lp);

    void findFeasibleStartPoint(); // phase 1

    void setLPState(const std::vector<double> &lpState);

    // TEST STUFF
    int countFractionalPivCols();
    int infeasibilityCount();

protected:
    void processProposal(const ProposalPivot &proposal);
    ProposalPivot proposePivot();
    int proposeColumn();


    void toCanonicalState();
    std::vector<int> auxiliaries();

    void updateLPSolution(const ProposalPivot &pivot);
    void revertLPSolution(const ProposalPivot &pivot);

    bool solutionIsPrimaryFeasible();

    double infeasibility();
};


#endif //GLPKTEST_SIMPLEXMCMC_H
