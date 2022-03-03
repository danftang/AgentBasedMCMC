//
// Created by daniel on 13/08/2021.
//

#ifndef GLPKTEST_EXPECTATIONPROBLEM_H
#define GLPKTEST_EXPECTATIONPROBLEM_H

// Represents the problem of finding a (vector) expectation value of a function, f,
// over a ConvexPMF. f should take a exactEndState to the ConvexPMF and return a vector
// which is an instance of the value we wish to calculate the expectation of.
//
// In order to calculate the expectation value, send an instancce of this class to a solver.
class ExpectationProblem {
    ConvexPMF pmf;
    std::function<std::vector<double>(const std::vector<double> &)> expectationFunction;

    ExpectationProblem(ConvexPMF pmf, std::function<std::vector<double>(const std::vector<double> &)> expectationFunc)
    :pmf(std::move(pmf)),
    expectationFunction(std::move(expectationFunc)) {

    }
};


#endif //GLPKTEST_EXPECTATIONPROBLEM_H
