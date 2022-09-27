#include "Experiments.h"
#include "FiguresForPaper.h"

int main(int argc, char *argv[]) {
//    auto startTime = std::chrono::steady_clock::now();

//    std::vector<int>::con it;

    SparseVec<double> vec;
    const SparseVec<double> &cvec = vec;

    vec.insert(3,5.678);
    vec.insert(7,6.789);

    for(const auto &entry: cvec) {
        std::cout << entry.index() << " " << entry.value() << std::endl;
    }

    for(double &value : vec.values) value *= 2.0;


    for(auto it = cvec.rbegin(); it != cvec.rend(); ++it) {
        std::cout << it->index() << " " << it->value() << std::endl;
    }

    /////////////// Generate data for figures 3 and 4 and statistics for table 2
//    FiguresForPaper<8,4>::generateStandardProblemFile(10.0);
//    FiguresForPaper<8,4>::generateStats(1000);

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
//    Experiments::CatMouseSingleObservation();
//    Experiments::CatMouseAssimilation();
//    Experiments::PredPreySingleObservation();
        // Trajectory: 8.16 x 59.56 = 486, RMS= 0.023, t= 58.8s.
        // Extended: 2.3 x 16.9 = 39.09, RMS = 0.030, sample time = 7.54s, 59.2% infeasible
        // Extended, separated phi factorial: 2.5 x 10.9 = 27.2, RMS=0.037, sample time = 7.6s, 62.9% infeasible
        // Extended, fully factorized 4.6 x 6.7 = 30.5, RMS=0.027, sample time = 2.8s, 53% infeasible
        // Fully extended, fully factorised 4.6x6.7 = 30.5, RMS=0.037, sample time = 2.8s(unplugged), 55.5% infeasible
//    auto endTime = std::chrono::steady_clock::now();
//    std::cout << "Exec time = " << endTime - startTime << std::endl;
    return 0;
}

