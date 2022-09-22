#include "Experiments.h"
#include "agents/BinomialAgent.h"
#include "subscript_operator_traits.h"
#include "BernoulliStartState.h"
#include "FiguresForPaper.h"


int main(int argc, char *argv[]) {
//    auto startTime = std::chrono::steady_clock::now();

    /////////////// Generate data for figures 3 and 4 and statistics for table 2
//    FiguresForPaper<32,16>::generateStandardProblemFile(10.0);
//    FiguresForPaper<32,16>::generateStats(10000000);
    /////////////// Plot figures 3 and 4 and print statistics for table 2
//    FiguresForPaper<32,16>::plotStats(false); // set to true to allow printing from plots.

    /////////////// Get timing for a sample on a single core
//    FiguresForPaper<32,16>::sampleTiming(1000000);

    ////////////// Test MCMC against small, tractable examples

    // TODO: try factorising extended model state

//    Experiments::Animation();

//    Experiments::BinomialAgentSingleObservation();
    Experiments::CatMouseSingleObservation();
//    Experiments::CatMouseAssimilation();
//    Experiments::PredPreySingleObservation();
        // Trajectory: 8.16 x 59.56 = 486, RMS= 0.023, t= 58.8s.
        // Extended: 2.3 x 16.9 = 39.09, RMS = 0.030, sample time = 7.54s, 59.2% infeasible
        // Extended, separated phi factorial: 2.5x10.9=27.2, RMS=0.037, sample time = 7.6s, 62.9% infeasible

//    auto endTime = std::chrono::steady_clock::now();
//    std::cout << "Exec time = " << endTime - startTime << std::endl;
    return 0;
}

