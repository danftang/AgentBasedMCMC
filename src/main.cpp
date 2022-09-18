#include "Experiments.h"
#include "agents/BinomialAgent.h"
#include "subscript_operator_traits.h"
#include "BernoulliStartState.h"
#include "FiguresForPaper.h"


int main(int argc, char *argv[]) {
    auto startTime = std::chrono::steady_clock::now();

    /////////////// Generate data for figures 3 and 4 and statistics for table 2
//    FiguresForPaper<32,16>::generateStandardProblemFile(10.0);
//    FiguresForPaper<32,16>::generateStats(10000000);
    /////////////// Plot figures 3 and 4 and print statistics for table 2
//    FiguresForPaper<32,16>::plotStats(false); // set to true to allow printing from plots.

    /////////////// Get timing for a sample on a single core
//    FiguresForPaper<32,16>::sampleTiming(1000000);

    ////////////// Test MCMC against small, tractable examples

    // TODO: try factorising extended model state

//    Experiments::BinomialAgentSingleObservation();
//    Experiments::CatMouseSingleObservation();
//    Experiments::CatMouseAssimilation();
    Experiments::PredPreySingleObservation(); // dep col = 8.16 dep row = 59.56, 29.7s

    auto endTime = std::chrono::steady_clock::now();
    std::cout << "Exec time = " << endTime - startTime << std::endl;
    return 0;
}

