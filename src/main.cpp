#include <iostream>
#include <future>
#include <optional>
#include <vector>

#include "Experiments.h"
#include "FiguresForPaper.h"
#include "include/asyncvector.h"
#include "agents/ReducedPredPreyAgent.h"
#include "agents/CatMouseAgent.h"
#include "TableauEntropyMaximiser.h"

int main(int argc, char *argv[]) {

    // Forward basis (8x4)
    // Entropy = 6.05847
    // mean basis size = 4.4
    // mean dependency col size = 5.22
    // mean dependency row size = 7.533
    // mean updates per transition = 39.35
    // infeasible = 69%
    // Sample inefficiency 4972, 1488
    // time per feasible sample 41.6uS
    // time per effective sample 207mS
    //
    // Backward basis (8x4)
    // Entropy = 5.97
    // mean basis size = 7.33
    // mean dependency col size = 11.03
    // mean dependency row size = 15.9
    // mean updates per transition = 175.3
    // infeasible = 62%
    // Sample inefficiency 123, 50.8
    // time per feasible sample 137uS
    // time per effective sample 16.8mS
    //
    // Norm-minimised backward-basis (8x4)
    // Entropy = 5.38
    // mean basis size = 7.33
    // mean dependency col size = 11.03
    // mean dependency row size = 16.0
    // mean updates per transition = 175.3
    // infeasible = 62%
    // Sample inefficiency = 91.9, 47.0
    // time per feasible sample 115.6uS
    // time per effective sample 10.6mS
    //

    // 8 x 4
    // kappa    %age infeasible     time/effetive sample
    // 7.25     53%                 113ms
    // 6.75     74%                 103ms
    // 6.5      81%                 75.9ms / 117.72
    // 6.25     88%                 94.5ms

    // No decay on start state or likelihood at negative vals: kappa=6.25, infeasible = 90%, time/eff-sample = 77ms
    // double decay on start state and likelihood: kappa=6.25, infeasible = 88%, time/eff-sample = 69.5ms
    // normal decay on start state, no decay on likelihood, marginal with step down on infeasible proir, kappa=6.5, infeasible=80%, time/eff-sample = 55.1ms

//    FiguresForPaper<32,16>::generateStandardPredPreyPosteriorFile(6.0);

    FiguresForPaper<8,4>::generateStandardPredPreyPosteriorFile(12.0);
    FiguresForPaper<8,4>::generateStats(200000);
    FiguresForPaper<8,4>::plotStats(false); // set to true to allow printing from plots.

    // 16 x 4
    // kappa    %age infeasible     time/effetive sample
    // 8.0      75%                 3035ms

//    FiguresForPaper<16,4>::generateStandardProblemFile(7.5);
//    FiguresForPaper<16,4>::generateStats(1000000);
//    FiguresForPaper<16,4>::plotStats(false); // set to true to allow printing from plots.

    /////////////// Generate data for figures 3 and 4 and statistics for table 2
    //    FiguresForPaper<32,16>::generateStandardProblemFile(10.0);
//    FiguresForPaper<32,16>::generateStats(10000000);
    /////////////// Plot figures 3 and 4 and print statistics for table 2
//    FiguresForPaper<32,16>::plotStats(false); // set to true to allow printing from plots.

    /////////////// Get timing for a sample on a single core
//    FiguresForPaper<32,16>::sampleTiming(1000000);

    ////////////// Test MCMC against small, tractable examples

//    Experiments::Animation();

// TODO: function entropies should always be between 0 and 1

//    Experiments::BinomialAgentSingleObservation();
//    Experiments::CatMouseSingleObservation();
//    Experiments::CatMouseAssimilation();
//    Experiments::PredPreySingleObservation();
//    Experiments::PredPreyPriorTest();
//    Experiments::ReducedPredPreyPriorTest();

        // Trajectory: 8.16 x 59.56 = 486, RMS= 0.023, t= 58.8s.
        // Extended: 2.3 x 16.9 = 39.09, RMS = 0.030, sample time = 7.54s, 59.2% infeasible
        // Extended, separated phi factorial: 2.5 x 10.9 = 27.2, RMS=0.037, sample time = 7.6s, 62.9% infeasible
        // Extended, fully factorized 4.6 x 6.7 = 30.5, RMS=0.027, sample time = 2.8s, 53% infeasible
        // Fully extended, fully factorised 4.6x6.7 = 30.5, RMS=0.037, sample time = 2.8s(unplugged), 55.5% infeasible
//    auto endTime = std::chrono::steady_clock::now();
//    std::cout << "Exec time = " << endTime - startTime << std::endl;
    return 0;
}

