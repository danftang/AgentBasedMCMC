// A distribution defined inside a convex polyhedron
// L <= CX <= H
// with a probability of the form
//
// P(X) = P_i(X) W(X)
//
// Where P_i(X) is factorized into functions of the auxiliary
// variables. i.e. expressed in the form
//
// P_i(X) = \prod_i F_i(C_iX)
//
// and W(X) is an importance weight function that can be updated
// by perturbations.
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_WEIGHTEDFACTOREDCONVEXDISTRIBUTION_H
#define ABMCMC_WEIGHTEDFACTOREDCONVEXDISTRIBUTION_H

#include "FactoredConvexDistribution.h"

template<class T, class IMPORTANCEFUNC>
class WeightedFactoredConvexDistribution: public FactoredConvexDistribution<T> {
public:
    IMPORTANCEFUNC importanceWeight;
};


#endif //ABMCMC_WEIGHTEDFACTOREDCONVEXDISTRIBUTION_H
