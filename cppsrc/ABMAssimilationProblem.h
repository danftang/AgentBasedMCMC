//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_ABMASSIMILATIONPROBLEM_H
#define GLPKTEST_ABMASSIMILATIONPROBLEM_H

template<typename AGENT, typename STARTSTATEPMF, typename LIKELIHOODPMF>
class ABMAssimilationProblem {
    // traits
    typedef AGENT Agent;
    typedef STARTSTATEPMF StartState;
    typedef LIKELIHOODPMF Likelihood;

    StartState startStatePMF;
    Likelihood likelihoodPMF;
    int nTimesteps;

    ABMAssimilationProblem(STARTSTATEPMF startState, LIKELIHOODPMF likelihood, int nTimesteps):
    startStatePMF(std::move(startState)),
    likelihoodPMF(std::move(likelihood)),
    nTimesteps(nTimesteps) {

    }
};


#endif //GLPKTEST_ABMASSIMILATIONPROBLEM_H
