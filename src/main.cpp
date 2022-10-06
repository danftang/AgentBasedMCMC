#include <iostream>
#include <future>
#include <optional>
#include <vector>

#include "Experiments.h"
#include "FiguresForPaper.h"

int main(int argc, char *argv[]) {


//    auto startTime = std::chrono::steady_clock::now();

//    std::cout
//            << exp(PredPreyAgentBase::lpPredBirthGivenPrey) << " "
//            << exp(PredPreyAgentBase::lpPredDeathGivenPrey) << " "
//            << exp(PredPreyAgentBase::lpPredBirthGivenNoPrey) << " "
//            << exp(PredPreyAgentBase::lpPredDeathGivenNoPrey) << " "
//            << exp(PredPreyAgentBase::lpPreyBirthGivenPred) << " "
//            << exp(PredPreyAgentBase::lpPreyDeathGivenPred) << " "
//            << exp(PredPreyAgentBase::lpPreyBirthGivenNoPred) << " "
//            << exp(PredPreyAgentBase::lpPreyDeathGivenNoPred) << std::endl;

//    std::cout
//            << exp(PredPreyAgentBase::lpPredDeath) << " "
//            << exp(PredPreyAgentBase::lpPredBirth) << " "
//            << exp(PredPreyAgentBase::lpPreyDeath) << " "
//            << exp(PredPreyAgentBase::lpPreyBirth) << std::endl;

    // 8 x 4
    // kappa    %age infeasible     time/effetive sample
    // 7.25     53%                 113ms
    // 6.75     74%                 103ms
    // 6.5      81%                 75.9ms / 117.72
    // 6.25     88%                 94.5ms

    // No decay on start state or likelihood at negative vals: kappa=6.25, infeasible = 90%, time/eff-sample = 77ms
    // double decay on start state and likelihood: kappa=6.25, infeasible = 88%, time/eff-sample = 69.5ms
    // normal decay on start state, no decay on likelihood, marginal with step down on infeasible proir, kappa=6.5, infeasible=80%, time/eff-sample = 55.1ms

//    FiguresForPaper<8,4>::generateStandardPredPreyPosteriorFile(6.0);
//    FiguresForPaper<8,4>::generateStats(100000,1);
//    FiguresForPaper<8,4>::plotStats(false); // set to true to allow printing from plots.


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

//    Experiments::BinomialAgentSingleObservation();
//    Experiments::CatMouseSingleObservation();
//    Experiments::CatMouseAssimilation();
    Experiments::PredPreySingleObservation();
        // Trajectory: 8.16 x 59.56 = 486, RMS= 0.023, t= 58.8s.
        // Extended: 2.3 x 16.9 = 39.09, RMS = 0.030, sample time = 7.54s, 59.2% infeasible
        // Extended, separated phi factorial: 2.5 x 10.9 = 27.2, RMS=0.037, sample time = 7.6s, 62.9% infeasible
        // Extended, fully factorized 4.6 x 6.7 = 30.5, RMS=0.027, sample time = 2.8s, 53% infeasible
        // Fully extended, fully factorised 4.6x6.7 = 30.5, RMS=0.037, sample time = 2.8s(unplugged), 55.5% infeasible
//    auto endTime = std::chrono::steady_clock::now();
//    std::cout << "Exec time = " << endTime - startTime << std::endl;
    return 0;
}

