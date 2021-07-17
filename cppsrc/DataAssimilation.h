//
// Created by daniel on 15/07/2021.
//

#ifndef GLPKTEST_DATAASSIMILATION_H
#define GLPKTEST_DATAASSIMILATION_H

#include <vector>
#include "ABMProblem.h"
#include "Trajectory.h"
#include "PoissonState.h"

template <typename AGENT>
class DataAssimilation {
public:

    // Calculates the information gained by assimilation over nTimesteps
    static double calculateInformationGain(
            int nWindows,
            int nTimestepsPerWindow,
            const PoissonState<AGENT> &startStatePrior,
            int nSamplesPerWindow
            ) {

    }


    static PoissonState<AGENT> assimilateWindow(
            int nTimeseps,
            std::vector<Observation<AGENT>> observations,
            const PoissonState<AGENT> &prior,
            const Trajectory<AGENT> &initialTrajectory,
            int nSamples) {

        ABMProblem<AGENT> abm(nTimeseps, observations, [&](const Trajectory<AGENT> &trajectory) {
            return prior.logProb(trajectory(0));
        });
        SimplexMCMC mcmc(abm, abm.logProbFunc());
        mcmc.setLPState(initialTrajectory);

        // test
        std::vector<int> occupationHistogram(10,0.0);

        PoissonState<AGENT> finalState;
        for(int n=0; n<nSamples; ++n) {
            mcmc.nextSample();
            if(nSamples<1000 || n%1000 == 1) {
//                std::cout << "Sample " << n << std::endl;
                assert(abm.isValidSolution(mcmc.X()));
//            std::cout << "Sample " << n << " : " << glp::SparseVec(mcmc.X()) << std::endl;
            }
            const Trajectory<AGENT> &trajectory = reinterpret_cast<const Trajectory<AGENT> &>(mcmc.X());
            ModelState<AGENT> endStateSample = trajectory(abm.nTimesteps);
            occupationHistogram[endStateSample[AGENT(152)]] += 1;
            finalState += endStateSample;
        }
        std::cout << "Occupation histogram " << occupationHistogram << std::endl;
        std::cout << "infeasible/feasible: " << mcmc.infeasibleStatistics.nSamples *100.0/mcmc.feasibleStatistics.nSamples << "%" << std::endl;
        std::cout << "Feasible sample statistics:" << std::endl << mcmc.feasibleStatistics << std::endl;
        std::cout << "Infeasible sample statistics:" << std::endl << mcmc.infeasibleStatistics << std::endl;
        return finalState;
    }


    // Calculates the posiion states at the end of each window, given no observations but
    // including the prior over the start state. This is needed in order to calculate the
    // information gain due to assimilation.
    static std::vector<PoissonState<AGENT>> prior(int nWindows, int nTimestepsPerWindow, int nSamples, const PoissonState<AGENT> &startState) {
        int nTimesteps = nWindows*nTimestepsPerWindow;
        std::vector<PoissonState<AGENT>> endOfWindowStates(nWindows);
        for(int s=0; s<nSamples; ++s) {
            Trajectory<AGENT> trajectory = Trajectory<AGENT>::run(startState, nTimesteps);
            for(int window=0; window<nWindows; ++window) {
                endOfWindowStates[window] += trajectory[(window+1)*nTimestepsPerWindow];
            }
        }
        return endOfWindowStates;
    }

    static double  informationGain(const ModelState<AGENT> &realState, const PoissonState<AGENT> &prior, const PoissonState<AGENT> &posterior) {
        return (posterior.logProb(realState) - prior.logProb(realState))/log(2);
    }


};


#endif //GLPKTEST_DATAASSIMILATION_H
