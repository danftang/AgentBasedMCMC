#include "Experiments.h"
#include "agents/CatMouseAgent.h"
// #include "FiguresForPaper.h"
#include "subscript_operator_traits.h"


int main(int argc, char *argv[]) {

    /////////////// Generate data for figures 3 and 4 and statistics for table 2
//    FiguresForPaper<32,16>::generateStandardProblemFile(10.0);
//    FiguresForPaper<32,16>::generateStats(10000000);
    /////////////// Plot figures 3 and 4 and print statistics for table 2
//    FiguresForPaper<32,16>::plotStats(false); // set to true to allow printing from plots.

    /////////////// Get timing for a sample on a single core
//    FiguresForPaper<32,16>::sampleTiming(1000000);

    ////////////// Test MCMC against small, tractable examples

    Experiments::CatMouseAssimilation();

//    Experiments::PredPreySingleObservation();
//    Experiments::CatMouseSingleObservation();
//    Experiments::BinomialAgentSingleObservation();

    return 0;
}

