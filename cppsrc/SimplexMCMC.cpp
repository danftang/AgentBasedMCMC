//
// Created by daniel on 15/04/2021.
//

#include "SimplexMCMC.h"
#include <algorithm>
#include <cmath>

SimplexMCMC::SimplexMCMC(glp::Problem &prob, const glp::SparseVec &initialSample) :
        Simplex(prob),
        colLastNonZero(n,-1),
        rowLatestCompletionPivot(m, -1),
        latestCompletionBegin(m, -1),
        lnRowPivotCount(m,0.0) {

    for(int i = 1; i<=m; ++i) {
        currentBasis[head[i]] = i;
    }
}


void SimplexMCMC::randomPivot() {

}

// Calculate the degeneracy probability of the current state
// from a cold start (i.e. nothing already calculated)
double SimplexMCMC::lnDegeneracyProb() {
    double row[n+1];
    int nextVarBegin;
    int nextVarEnd = n+m+1;
    int i,j,k,last,lastj;
    int c = 0;
    double lnP = 0.0;

    auto basisVar = currentBasis.rbegin();
    while(basisVar != currentBasis.rend()) {
        i = basisVar->second;
        tableauRow(i, row);
        ++basisVar;
        nextVarBegin = basisVar->first + 1;
        last = 0;
        for(j = 1; j <=n; ++j) {
            if(row[j] != 0.0) {
                k = head[m+j];
                if(colLastNonZero[j] == -1) {
                    colLastNonZero[j] = i;
                    if(k >= nextVarBegin && k < nextVarEnd) ++c;
                }
                if(k>last) {
                    last = k;
                    lastj = j;
                }
            }
        }
        rowLatestCompletionPivot[i] = last;
        lnP += lnRowPivotCount[i] = std::log(c);
        if(last < nextVarEnd) nextVarEnd = last;
        latestCompletionBegin[i] = nextVarEnd;
        pivot(i,lastj,false);
    }
    return lnP;
}

