#include "Experiments.h"
#include "agents/BinomialAgent.h"
// #include "FiguresForPaper.h"
#include "subscript_operator_traits.h"
#include "BernoulliStartState.h"


int main(int argc, char *argv[]) {

//    BernoulliStartState<BinomialAgent<3>> startState({1.0, 0.1, 0.0});
//    std::cout << startState << std::endl;
//    auto sampler = startState.sampler();
//    for(int s=0; s<100; ++s) std::cout << sampler() << std::endl;

    /////////////// Generate data for figures 3 and 4 and statistics for table 2
//    FiguresForPaper<32,16>::generateStandardProblemFile(10.0);
//    FiguresForPaper<32,16>::generateStats(10000000);
    /////////////// Plot figures 3 and 4 and print statistics for table 2
//    FiguresForPaper<32,16>::plotStats(false); // set to true to allow printing from plots.

    /////////////// Get timing for a sample on a single core
//    FiguresForPaper<32,16>::sampleTiming(1000000);

    ////////////// Test MCMC against small, tractable examples

//    Experiments::BinomialAgentSingleObservation();
//    Experiments::CatMouseSingleObservation();
//    Experiments::CatMouseAssimilation();
//    Experiments::PredPreySingleObservation();

    return 0;
}

