#include <iostream>
#include "Experiments.h"
#include "FiguresForPaper.h"

//void testActFermionicDistribution() {
//    std::vector<double> p = {0.8, 0.15, 0.05};
//    ActFermionicDistribution myDist(p);
//
//    std::vector<int> actCounts(3,0);
//    for(int s = 0; s<1000000; ++s) {
//        std::vector<bool> acts = myDist.sampleUnordered(2);
//        if(acts[0]) ++actCounts[0];
//        if(acts[1]) ++actCounts[1];
//        if(acts[2]) ++actCounts[2];
//    }
//    std::cout << actCounts << std::endl;
//
//    std::vector<int> actCounts2(3,0);
//    std::discrete_distribution<int> dist = std::discrete_distribution<int>(p.begin(), p.end());
//    for(int s = 0; s<1000000; ) {
//        int act1 = dist(Random::gen);
//        int act2 = dist(Random::gen);
//        if(act1 != act2) {
//            ++actCounts2[act1];
//            ++actCounts2[act2];
//            ++s;
//        }
//    }
//    std::cout << actCounts2 << std::endl;
//    assert(abs(actCounts[0]-actCounts2[0]) < 1000);
//    assert(abs(actCounts[1]-actCounts2[1]) < 1000);
//    assert(abs(actCounts[2]-actCounts2[2]) < 1000);
//}


int main(int argc, char *argv[]) {

//    Experiments::CatMouseAssimilation();

//    Experiments::PredPreySingleObservation();
//    Experiments::CatMouseSingleObservation();
//    Experiments::BinomialAgentSingleObservation();

//    FiguresForPaper::generateStandardProblemFile<8>(4, 6.0);
//    FiguresForPaper::singleThreadStats<8>(4,100000);
//    FiguresForPaper::generateStatsAndPlot<8>(4, 100000);

//    FiguresForPaper::generateStandardProblemFile<8>(8, 6.5);
//    FiguresForPaper::singleThreadStats<8>(8,100000);
//    FiguresForPaper::generateStatsAndPlot<8>(8, 500000);
//    FiguresForPaper::plotStats<8>(8);

//    FiguresForPaper::generateStandardProblemFile<8>(16, 7.5);
//    FiguresForPaper::singleThreadStats<8>(16,100000);
//    FiguresForPaper::generateStatsAndPlot<8>(16, 500000);
//    FiguresForPaper::plotStats<8>(16);


//    FiguresForPaper::generateStandardProblemFile<10>(5, 6.5);
//    FiguresForPaper::singleThreadStats<10>(5,100000);
//    FiguresForPaper::generateStatsAndPlot<10>(5, 200000);
//    FiguresForPaper::plotStats<10>(5);

//    FiguresForPaper::generateStandardProblemFile<16>(4, 6.0);
//    FiguresForPaper::singleThreadStats<16>(4,100000);
//    FiguresForPaper::generateStatsAndPlot<16>(4,100000);
//    FiguresForPaper::plotStats<16>(4);

//    FiguresForPaper::generateStandardProblemFile<16>(8, 7.75);
//    FiguresForPaper::singleThreadStats<16>(8,100000);
//    FiguresForPaper::generateStatsAndPlot<16>(8, 50000);
//    FiguresForPaper::plotStats<16>(8);

//    FiguresForPaper::generateStandardProblemFile<16>(16, 9.0);
//    FiguresForPaper::singleThreadStats<16>(16,100000);
//    FiguresForPaper::generateStatsAndPlot<16>(16, 4000000);
//    FiguresForPaper::plotStats<16>(16);


//    FiguresForPaper::generateStandardProblemFile<32>(4, 9.5);
//    FiguresForPaper::singleThreadStats<32>(4,100000);
//    FiguresForPaper::generateStatsAndPlot<32>(4,3000000);

//    FiguresForPaper::generateStandardProblemFile<32>(8, 10.0);
//    FiguresForPaper::singleThreadStats<32>(8,100000);
//    FiguresForPaper::generateStatsAndPlot<32>(8,5000000);

//    FiguresForPaper::generateStandardProblemFile<32>(16, 10.5);
//    FiguresForPaper::singleThreadStats<32>(16,100000);
//    FiguresForPaper::generateStatsAndPlot<32>(16, 12000000);
    FiguresForPaper::plotStats<32>(16,true);

//    FiguresForPaper::plotProblemEndState<32>(16);

//    FiguresForPaper::generateAllProblemFiles();


//    Experiments::minimalBasis();

//    Experiments::animatedPredPreyDemo();

//    FiguresForPaper::generateStatsAndPlot<8>(8);
//    Experiments::PredPreyAssimilation();


//    FiguresForPaper::plotStats<8>(2);
//    FiguresForPaper::plotStats<8>(4);
//    FiguresForPaper::plotStats<8>(6);
//    FiguresForPaper::plotStats<8>(8);
//   FiguresForPaper::plotStats<16>(8);
//    FiguresForPaper::plotStats<32>(8);
//    FiguresForPaper::plotStats<32>(8);
//    FiguresForPaper::plotStats<32>(16);


//    FiguresForPaper::generateStandardProblemFile();
//    Experiments::PredPreyConvergence();

//    Experiments::animatedPredPreyDemo();
//        Experiments::DataflowDemo();
//    Experiments::CatMouseMultiObservation();
//    Experiments::PredPreyAssimilation();

//    Experiments::FermionicIntegrality();

//    PoissonState<PredPreyAgent> startState;
//    TrajectoryPriorDistribution<PredPreyAgent> myPrior(startState, 8);

//    UnitTests tests;
//    tests.testPriorSampler();
//    tests.testRejectionSampler();
//    tests.testABMPrior();
//    tests.testSimplexSampler();
//    tests.testExactSolver();

//        std::cout << Experiments::informationIncrease(
//            8,
//            2,
//            1,
//            0.16,
//            0.32,
//            0.05,
//            0.9,
//            100000,
//            1000) << std::endl;


//    Experiments::RandomWalk();
//    Experiments::GnuplotTest();

    return 0;
}

