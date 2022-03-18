// A distribution defined inside a convex polyhedron
// L <= CX <= H
// with a probability that can be factorized into functions
// of the auxiliary variables. i.e. expressed in the form
//
// P(X) = \prod_i F_i(C_iX)
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_FACTOREDCONVEXDISTRIBUTION_H
#define ABMCMC_FACTOREDCONVEXDISTRIBUTION_H

#include <functional>
#include "ConvexPolyhedron.h"

template<class T>
class FactoredConvexDistribution {
public:
    ConvexPolyhedron<T>                     constraints;
    std::vector<std::function<double(T)>>   factors;

    void addFactor(Constraint<T> constraint, std::function<double(T)> factor) {
        constraints.push_back(std::move(constraint));
        factors.push_back(std::move(factor));
    }
};


#endif //ABMCMC_FACTOREDCONVEXDISTRIBUTION_H
