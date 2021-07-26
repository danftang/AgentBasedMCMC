//
// Created by daniel on 20/07/2021.
//

#ifndef GLPKTEST_ASSIMILATIONWINDOW_H
#define GLPKTEST_ASSIMILATIONWINDOW_H

#include <vector>
#include "debug.h"

template<typename AGENT>
class AssimilationWindow {
public:
    Trajectory<AGENT>               realTrajectory;
    std::vector<Observation<AGENT>> observations;
    PoissonState<AGENT>             analysis;


    AssimilationWindow(const Trajectory<AGENT> &realTrajectory,
                       const PoissonState<AGENT> &priorStartState,
                       std::vector<Observation<AGENT>> observations,
                       int nSamples,
                       int burnInSamples) :
            realTrajectory(realTrajectory),
            observations(observations) {
        doAnalysis(priorStartState, nSamples, burnInSamples);
    }


    void doAnalysis(const PoissonState<AGENT> &priorStartState, int nSamples, int burnInSamples) {
        ABMProblem<AGENT> abm(realTrajectory.nTimesteps(), observations, [&](const Trajectory<AGENT> &trajectory) {
            return priorStartState.logProb(trajectory(0));
        });
        SimplexMCMC mcmc(abm, abm.logProbFunc());

        Trajectory<AGENT> firstGuessTrajectory(realTrajectory.nTimesteps(), priorStartState.sample());
        mcmc.setLPState(firstGuessTrajectory);
        mcmc.findFeasibleStartPoint();
        assert(mcmc.abmSanityChecks());

        for(int burnIn=0; burnIn < burnInSamples; ++burnIn) { mcmc.nextSample(); }
        for(int n=0; n<nSamples; ++n) {
            mcmc.nextSample();
            debug(if(nSamples<1000 || n%1000 == 1) {
                assert(abm.isValidSolution(mcmc.X()));
                std::cout << "Got sample " << mcmc.X() << std::endl;
            });
            const Trajectory<AGENT> &trajectory = reinterpret_cast<const Trajectory<AGENT> &>(mcmc.X());
            analysis += trajectory.endState();
        }
        debug(std::cout
                      << "infeasible/feasible: " << mcmc.infeasibleStatistics.nSamples *100.0/mcmc.feasibleStatistics.nSamples << "%" << std::endl
                      << "Feasible sample statistics:" << std::endl << mcmc.feasibleStatistics << std::endl
                      << "Infeasible sample statistics:" << std::endl << mcmc.infeasibleStatistics << std::endl;
        );
    }

    friend Gnuplot &operator<<(Gnuplot &gp, const AssimilationWindow<AGENT> &window) {
        typedef std::tuple<double, double, double, double, double> HeatRecord;
        std::vector<std::vector<HeatRecord>> heatData;
        std::vector<std::tuple<double, double, double>> pointData;

        ModelState<AGENT> realState = window.realTrajectory.endState();

        for (int x = 0; x < PredPreyAgent::GRIDSIZE; ++x) {
            for (int y = 0; y < PredPreyAgent::GRIDSIZE; ++y) {
                int colour = 2 * (realState[PredPreyAgent(x, y, PredPreyAgent::PREDATOR)] > 0.0)
                             + (realState[PredPreyAgent(x, y, PredPreyAgent::PREY)] > 0.0);
                if (colour != 0)
                    pointData.emplace_back(x, y, colour);
            }
        }

        for (int x = 0; x < PredPreyAgent::GRIDSIZE; ++x) {
            std::vector<HeatRecord> &record = heatData.emplace_back();
            for (int y = 0; y < PredPreyAgent::GRIDSIZE; ++y) {
                double lPrey = window.analysis.lambda(PredPreyAgent(x, y, PredPreyAgent::PREY));
                double lPred = window.analysis.lambda(PredPreyAgent(x, y, PredPreyAgent::PREDATOR));
                record.emplace_back(x, y, std::min(lPrey, 1.0) * 200.0, 0.0, std::min(lPred, 1.0) * 200.0);
            }
        }

        gp << "set linetype 1 lc 'red'\n";
        gp << "set linetype 2 lc 'blue'\n";
        gp << "set linetype 3 lc 'magenta'\n";
        gp << "plot [-0.5:" << PredPreyAgent::GRIDSIZE - 0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE - 0.5 << "] ";
        gp << "'-' with rgbimage notitle, ";
        gp << "'-' with points pointtype 5 pointsize 0.5 lc variable notitle\n";
        gp.send2d(heatData);
        gp.send1d(pointData);
        return gp;
    }
};


#endif //GLPKTEST_ASSIMILATIONWINDOW_H
