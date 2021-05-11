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
#include "BasisProbability.h"
#include "PivotPoint.h"

class SimplexMCMC: public glp::Simplex {
public:

    std::vector<int>    colLastNonZero;             // row-index of the last non-zero entry in this col
    std::vector<int>    rowLatestCompletionPivot; // k-col-index of the earliest col of the final basis that can pivot on this row
    std::vector<int>    latestCompletionBegin;   // k-col-index of the earliest col of the final completion
    std::vector<double> lnRowPivotCount;            // ln of number of possible choices of next in-sequence var of degenerate state
    std::set<int>       finalBasis;                 // the final degenerate state (ordered)
    std::map<int,int>   currentBasis;               // the current basis ordered map from k-index to m-index
    int nSamples = 0;
    int nRejections = 0;
    int fractionalRunLength = 0;

    BasisProbability probability;


    SimplexMCMC(glp::Problem &prob);

    void randomPivot();
    double lnDegeneracyProb();

    glp::SparseVec nextSample();
    void processProposal(PivotPoint proposal);
    double transitionProb(PivotPoint proposal);
    double reverseTransitionProb(PivotPoint proposal);

    PivotPoint proposePivot();

    std::vector<int> calcPivotRows(int j, const std::vector<double> &colVec);
};


#endif //GLPKTEST_SIMPLEXMCMC_H
