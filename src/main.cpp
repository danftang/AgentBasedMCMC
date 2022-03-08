#include <iostream>
#include <thread>
#include <future>
//#include "agents/CatMouseAgent.h"
//#include "PredPreyProblem.h"
//#include "Trajectory.h"
//#include "Experiments.h"
//#include "StlStream.h"
//#include "Random.h"
//#include "UnitTests.h"
//#include "ConvexPolyhedron.h"
//#include "ABMConstraints.h"
//#include "diagnostics/Dataflow.h"
//#include "diagnostics/MeanAndVariance.h"
//#include "Plotter.h"
//#include "FiguresForPaper.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include "agents/BinomialAgent.h"
#include "BernoulliModelState.h"
#include "AssimilationWindow.h"
#include "TableauNormMinimiser.h"
#include "SparseBasis.h"

void BinomialAgentAssimilation() {
    constexpr int nTimesteps = 2;
    BinomialAgent::GRIDSIZE = nTimesteps + 1;
    constexpr int nSamplesPerWindow = 1000000;
    constexpr int nBurninSamples = 5000;

    BernoulliModelState<BinomialAgent> startState([](BinomialAgent agent) {
        return agent.stateId==0?1.0:0.1;
    });

    AgentStateObservation<BinomialAgent> observation(State<BinomialAgent>(1, 0),1,0.9);

    AssimilationWindow<BinomialAgent> window(nTimesteps, startState, observation);
    std::cout << "Prior support is \n" << window.priorPMF.convexSupport << std::endl;
    std::cout << "Likelihood support is \n" << window.likelihoodPMF.convexSupport << std::endl;
    std::cout << "Posterior support is \n" << window.posterior.convexSupport << std::endl;

    TableauNormMinimiser tableau(window.posterior.convexSupport);

    std::cout << "Initial tableau" << std::endl << tableau << std::endl;
    tableau.findMinimalBasis();
    std::cout << "Reduced tableau" << std::endl << tableau << std::endl;
    SparseBasis<int> basis(tableau);
    std::cout << "Basis is" << std::endl << basis << std::endl;

//    MCMCSampler<Trajectory<BinomialAgent>> sampler(window.posterior, window.priorSampler());
//    std::cout << "simplex = \n" << sampler.simplex << std::endl;
//    std::cout << "Starting burn-in" << std::endl;
//    for(int s=0; s<nBurninSamples; ++s)  sampler();
//    std::cout << "Done burn-in" << std::endl;
//    ModelStateSampleStatistics<BinomialAgent> sampleStats(sampler, nSamplesPerWindow);

//    std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
//    std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
//    std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;
//    std::vector<double> analysis = sampleStats.means();
//    std::cout << "Analysis means = " << analysis << std::endl;

//    std::cout << "Starting exact solve" << std::endl;
//    ExactSolver<BinomialAgent> exactSolver(window.posterior);
//    std::cout << "Exact means = " << exactSolver.exactEndState << std::endl;
//    double mse = 0.0;
//    for(int i=0; i<analysis.size(); ++i) {
//        mse += pow(analysis[i] - exactSolver.exactEndState[i],2.0);
//    }
//    mse /= analysis.size();
//    std::cout << "Mean square error = " << mse << std::endl;
}


int main(int argc, char *argv[]) {

    boost::numeric::ublas::vector<int> X(4, 1);
    boost::numeric::ublas::compressed_vector<int> Y(32000000000);
    boost::numeric::ublas::mapped_vector<int> Z(4);
//    Z.insert_element(1, 1);
    std::cout << Z.size() << std::endl;
    std::cout << Y.size() << std::endl;
    Z[1] = 1;
    Y = Z;
    Y.resize(4);
    Y.reserve(4);
    Y.insert_element(3,2);
    auto it = Y.begin();
    for(auto it= Y.begin(); it != Y.end(); it++) {
        std::cout << it.index() << ", " << *it <<  std::endl;
    }

    X = X + Y;
    std::cout << Z[0] << ", " << Z[1] << std::endl;

    std::cout << X[0] << ", " << X[1] << std::endl;


    BinomialAgentAssimilation();


//    Experiments::minimalBasis();

//    Experiments::animatedPredPreyDemo();
//    FiguresForPaper::generateAllProblemFiles();

//    FiguresForPaper::generateStatsAndPlot<8>(8);
//    Experiments::PredPreyAssimilation();

//    FiguresForPaper::generateStats<8>(8);
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
//  Experiments::CatMouseSingleObservation();
//    Experiments::CatMouseAssimilation();
//    Experiments::CatMouseMultiObservation();
//    Experiments::PredPreySingleObservation();
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

