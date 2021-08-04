//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_ABMREJECTIONSAMPLER_H
#define GLPKTEST_ABMREJECTIONSAMPLER_H

template<typename PROBLEM>
class ABMRejectionSampler {
    const PROBLEM &problem;

    ABMRejectionSampler(const PROBLEM &prob): problem(prob) {}

    // N.B. Only use this when likelihood is reasonably high
    Trajectory<PROBLEM::Agent> nextSample() {
        double logLikelihood;
        Trajectory<PROBLEM::Agent> sample(0);
        do {
            sample = Trajectory<PROBLEM::Agent>(
                    problem.nTimesteps,
                    ModelState<PROBLEM::Agent>(problem.startState.nextSample())
                    );
            logLikelihood = problem.likelihood(sample);
        } while(Random::nextDouble() > exp(logLikelihood));
        return sample;
    }

};


#endif //GLPKTEST_ABMREJECTIONSAMPLER_H
