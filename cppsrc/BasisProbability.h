//
// Created by daniel on 11/05/2021.
//

#ifndef GLPKTEST_BASISPROBABILITY_H
#define GLPKTEST_BASISPROBABILITY_H

#include "glpkpp.h"

class SimplexMCMC;

class BasisProbability {
public:
    glp::SparseVec currentSample;
    double logDegeneracyProb;
    double logPX;
    double logFractionalPenalty;

    double logProb() {
        return logPX + logDegeneracyProb + logFractionalPenalty;
    }

    BasisProbability(SimplexMCMC &simplex) { }

};


#endif //GLPKTEST_BASISPROBABILITY_H
