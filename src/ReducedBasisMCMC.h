// MCMC on a linear manifold within a hyper-rectangle
// with a fixed basis
//
// Created by daniel on 03/03/2022.
//

#ifndef GLPKTEST_REDUCEDBASISMCMC_H
#define GLPKTEST_REDUCEDBASISMCMC_H


#include <vector>

class ReducedBasisMCMC {

    std::vector<double> &operator()(); // get next sample
};


#endif //GLPKTEST_REDUCEDBASISMCMC_H
