//
// Created by daniel on 15/04/2021.
//

#ifndef GLPKTEST_SIMPLEXMCMC_H
#define GLPKTEST_SIMPLEXMCMC_H

#include <vector>
#include <set>
#include <map>
#include "glpkpp.h"

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

    SimplexMCMC(glp::Problem &prob, const glp::SparseVec &initialSample);

    void randomPivot();
    double lnDegeneracyProb();

//
//
//    fun<R> expectation(nSamples: Int, initialExpectation: R, expectationAccumulator: (SparseVector<Fraction>, R) -> R): R {
//        var e = initialExpectation
//        var oldSample: SparseVector<Fraction>? = null
//        var rejections = 0
//        var lastTime = Instant.now().toEpochMilli()
//        for(s in 1..nSamples) {
//            val newSample = simplex.nextIntegerSample()
//            e = expectationAccumulator(newSample, e)
//            if(oldSample === newSample) ++rejections
//            oldSample = newSample
//            if(s.rem(100) == 0) {
//                val now = Instant.now().toEpochMilli()
//                println("Got to sample $s in ${(now-lastTime)/1000.0}s largest Numerator,Denominator ${largestNumeratorDenominator()}, Size ${simplex.M.rows.sumBy { it.nonZeroEntries.size }}, Degeneracy ${simplex.degeneracy()} logPiv = ${simplex.state.logProbOfPivotState} logPX = ${simplex.state.logPX}, logPDegeneracy = ${simplex.state.logDegeneracyProb} prob per event = ${simplex.state.logPX/newSample.nonZeroEntries.size}")
//                lastTime = now
////                println(simplex.fractionalLogP - simplex.state.logPX - ln(simplex.transitionProb(simplex.proposePivot())))
//            }
//        }
//        println("Rejection ratio = ${rejections.toDouble()/nSamples}")
//        return e
//    }

};


#endif //GLPKTEST_SIMPLEXMCMC_H
