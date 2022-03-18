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
//#include "BernoulliModelState.h"
//#include "AssimilationWindow.h"
//#include "TableauNormMinimiser.h"
//#include "SparseBasisSampler.h"
#include "BernoulliStartState.h"
#include "ObservationLikelihood.h"
#include "TrajectoryForecastDistribution.h"

void BinomialAgentAssimilation() {
    constexpr int nTimesteps = 2;
    typedef BinomialAgent<nTimesteps+1> AGENT;
    constexpr int nSamplesPerWindow = 1000000;
    constexpr int nBurninSamples = 5000;

    BernoulliStartState<AGENT> PStartState([](AGENT agent) { return agent.stateId==0?1.0:0.1; });
    std::cout << "Start state support is \n" << PStartState.constraints << std::endl;

//    ObservationLikelihood<BinomialAgent> PObsGivenTrajectory(AgentStateObservation(State<BinomialAgent>(1, 0),1,0.9));
//    TrajectoryForecastDistribution<BinomialAgent> PTrajectoryGivenStartState(nTimesteps);
//
//    // WeightedFactoredConvexDistribution<BinomialAgent,int>
//    auto posterior = PObsGivenTrajectory * PTrajectoryGivenStartState * PStartState;
//
//    std::cout << "Likelihood support is \n" << PObsGivenTrajectory << std::endl;
//    std::cout << "Posterior support is \n" << posterior << std::endl;
//
//    SparseBasisSampler<int> sampler(posterior);

    //...

//    TableauNormMinimiser tableau(posterior);
//    std::cout << "Initial tableau" << std::endl << tableau << std::endl;
//    tableau.findMinimalBasis();
//    std::cout << "Reduced tableau" << std::endl << tableau << std::endl;
//    SparseBasisSampler<int> basis(tableau, [](int i, int occupation) {
//        Event<BinomialAgent> event(i);
//        return event.agent().marginalTimestep(event.act());
//    });
//    std::cout << "Basis is" << std::endl << basis << std::endl;

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

