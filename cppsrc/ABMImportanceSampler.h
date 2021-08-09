//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_ABMIMPORTANCESAMPLER_H
#define GLPKTEST_ABMIMPORTANCESAMPLER_H

// PROBLEM should be an Assimilation Problem with traits and members as in AssimilationProblem
template<typename PROBLEM>
class ABMImportanceSampler {
public:
    const PROBLEM &problem;

    ABMImportanceSampler(const PROBLEM &prob): problem(prob) {}

    // N.B. Weight is log and unnormalised
    std::pair<Trajectory<PROBLEM::Agent>,double> nextSample() {
        double logWeight;
        Trajectory<PROBLEM::Agent> sample(0);
        do {
            sample = Trajectory<PROBLEM::Agent>(
                    problem.nTimesteps,
                    ModelState<PROBLEM::Agent>(problem.startState.nextSample())
                            );
            logWeight = problem.logLikelihood(sample);
        } while(logWeight == -INFINITY);
        return std::pair(sample, logWeight);
    }
};


#endif //GLPKTEST_ABMIMPORTANCESAMPLER_H
