#include <iostream>
#include "Experiments.h"
#include "FiguresForPaper.h"


int main(int argc, char *argv[]) {

//    Experiments::CatMouseAssimilation();

//    Experiments::PredPreySingleObservation();
//    Experiments::CatMouseSingleObservation();
//    Experiments::BinomialAgentSingleObservation();

//    FiguresForPaper::generateStandardProblemFile<8>(4, 5.75, 2.0);
//    FiguresForPaper::singleThreadStats<8>(4);
//    FiguresForPaper::generateStats<8>(4);
    FiguresForPaper::generateStatsAndPlot<8>(4);


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

