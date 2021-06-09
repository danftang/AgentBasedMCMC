//
// Created by daniel on 11/05/2021.
//

#include "BasisProbability.h"
#include "SimplexMCMC.h"

BasisProbability::BasisProbability(SimplexMCMC &simplex) {
    logDegeneracyProb = 0.0;//simplex.lnDegeneracyProb();
    logPX = simplex.lnProb();
    logFractionalPenalty = simplex.lnFractionalPenalty();
}
