//
// Created by daniel on 06/05/2021.
//

#include <vector>
#include "Experiments.h"
#include "agents/CatMouseAgent.h"
#include "ABMProblem.h"
#include "SimplexMCMC.h"
#include "agents/PredPreyAgent.h"
#include "PoissonState.h"
#include "DataAssimilation.h"
#include "debug.h"
#include "BinomialDistribution.h"
#include "TrajectoryPriorDistribution.h"
#include "SampleStatistics.h"
#include "TrajectorySampler.h"
#include "TrajectoryLikelihoodPMF.h"
#include "ExactSolver.h"
#include "agents/BinomialAgent.h"
#include "MCMCSolver.h"
#include "ABMPlotter.h"
#include "StlStream.h"

std::vector<double> Experiments::informationIncrease(int argc, char *argv[]) {
    if(argc != 9) {
        std::cout << "Wrong number of arguments. Should be <GRIDSIZE> <nTimestepsPerWindow> <nWindows> <pPredator> <pPrey> <pMakeObservation> <pObserveIfPresent> <nSamplesPerWindow> <nBurnInSamples>" << std::endl;
        return std::vector<double>();
    }
    glp_term_out(GLP_OFF); // turn off GLPK terminal output
    int gridsize = atoi(argv[1]);
    int windowSize = atoi(argv[2]);
    int nWindows = atoi(argv[3]);
    double pPredator = atof(argv[4]);         // Poisson prob of predator in each gridsquare at t=0
    double pPrey = atof(argv[5]);             // Poisson prob of prey in each gridsquare at t=0
    double pMakeObservation = atof(argv[6]);  // prob of making an observation of each gridsquare at each timestep
    double pObserveIfPresent = atof(argv[7]); // prob of detecting an agent given that it is present
    int nSamplesPerWindow = atoi(argv[8]);
    int nBurnInSamples = atoi(argv[9]);
    return informationIncrease(gridsize, windowSize, nWindows, pPredator, pPrey,
                               pMakeObservation, pObserveIfPresent, nSamplesPerWindow, nBurnInSamples);
}


std::vector<double> Experiments::informationIncrease(
        int gridsize,
        int windowSize,
        int nWindows,
        double pPredator,         // Poisson prob of predator in each gridsquare at t=0
        double pPrey,             // Poisson prob of prey in each gridsquare at t=0
        double pMakeObservation,  // prob of making an observation of each gridsquare at each timestep
        double pObserveIfPresent, // prob of detecting an agent given that it is present
        int nSamplesPerWindow,
        int nBurnInSamples
        ) {
    glp_term_out(GLP_OFF); // turn off GLPK terminal output
    PredPreyAgent::GRIDSIZE = gridsize;

    std::vector<boost::math::binomial_distribution<double>> startWeights(PredPreyAgent::domainSize());
    for(int agentId=0; agentId < PredPreyAgent::domainSize(); ++agentId) {
        startWeights[agentId] = boost::math::binomial_distribution<double>(
                1,(PredPreyAgent(agentId).type() == PredPreyAgent::PREDATOR?pPredator:pPrey)
                );
    }
    BinomialDistribution startStateDist(startWeights);


    AssimilationWindow<PredPreyAgent> window(windowSize, startStateDist, pMakeObservation, pObserveIfPresent);

    MCMCSolver<PredPreyAgent> solver(window);
    solver.solve(nBurnInSamples, nSamplesPerWindow);


    return std::vector<double>(); //assimilation.calculateInformationGain();
}


void Experiments::PredPreyAssimilation() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    PredPreyAgent::GRIDSIZE = 8;
    constexpr int windowSize = 2;
    constexpr int nWindows = 1;
    constexpr double pPredator = 0.08;//0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
    constexpr double pMakeObservation = 0.1;//0.04;    // prob of making an observation of each gridsquare at each timestep
    constexpr double pObserveIfPresent = 0.95;
    constexpr int nSamplesPerWindow = 1000000; //250000;
    constexpr int nBurninSamples = 1000;
//    constexpr int plotTimestep = nTimesteps-1;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////

//    std::vector<boost::math::binomial_distribution<double>> startWeights(PredPreyAgent::domainSize());
//    for(int agentId=0; agentId < PredPreyAgent::domainSize(); ++agentId) {
//        startWeights[agentId] = boost::math::binomial_distribution<double>(
//                1,(PredPreyAgent(agentId).type() == PredPreyAgent::PREDATOR?pPredator:pPrey)
//                );
//    }
//    BinomialDistribution startStateDist(startWeights);

    IntSampleStatistics startStateDist(
            PredPreyAgent::domainSize(),
            [](int agentId, int count) {
                double p = (PredPreyAgent(agentId).type() == PredPreyAgent::PREDATOR ? pPredator : pPrey);
                switch (count) {
                    case 0: return 1.0 - p;
                    case 1: return p;
                }
                return 0.0;
            });

    std::cout << "Initial state distribution = " << startStateDist << std::endl;


    TrajectorySampler<PredPreyAgent> priorSampler(nWindows * windowSize, startStateDist.sampler());
//    priorSampler.endState(100000);
    for(int s=0; s<100000; ++s) { // TODO: Work out why this causes much higher occupancy
//        startStateDist.nextSample();
        priorSampler.nextSample();
    }

//    std::cout << "start state sample: " << startStateDist.nextSample() << std::endl;

    DataAssimilation<PredPreyAgent> assimilation(startStateDist, pMakeObservation, pObserveIfPresent);

//    std::cout << "assimilation start state sample: " << assimilation.analysisSampler() << std::endl; // TODO: and why this fixes it!


    for(int w=0; w<nWindows; ++w) {
        const AssimilationWindow<PredPreyAgent> &window = assimilation.addWindow(windowSize, nBurninSamples, nSamplesPerWindow);

        std::cout << "Analysis = " << assimilation.analysis << std::endl;
        std::cout << "Means = " << assimilation.analysis.means() << std::endl;

        ABMPlotter<PredPreyAgent> gp;
        gp.plot(window.realTrajectory.endState(), assimilation.analysis.means());
//        BinomialDistribution prior = window.priorEndState(100000);
//        std::cout << "Window information gain = "
//        << informationGain(window.realTrajectory.endState(), prior, assimilation.analysis) << std::endl;
    }

    IntSampleStatistics priorEnd = priorSampler.endState(100000);
        std::cout << "Prior end distribution = " << priorEnd << std::endl;

    std::cout << "Total information gain = "
    << informationGain(
            assimilation.windows.back().realTrajectory.endState(),
            priorEnd,
            assimilation.analysis) << std::endl;


}


void Experiments::PredPreyExpt() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    PredPreyAgent::GRIDSIZE = 3;
    constexpr int windowSize = 2;
    constexpr int nWindows = 1;
    constexpr double pPredator = 0.08;//0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
    constexpr double pMakeObservation = 0.1;//0.04;    // prob of making an observation of each gridsquare at each timestep
    constexpr double pObserveIfPresent = 0.95;
    constexpr int nSamplesPerWindow = 500000; //250000;
    constexpr int nBurninSamples = 10000;
    //    constexpr int plotTimestep = nTimesteps-1;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////

    //    std::vector<boost::math::binomial_distribution<double>> startWeights(PredPreyAgent::domainSize());
    //    for(int agentId=0; agentId < PredPreyAgent::domainSize(); ++agentId) {
    //        startWeights[agentId] = boost::math::binomial_distribution<double>(
    //                1,(PredPreyAgent(agentId).type() == PredPreyAgent::PREDATOR?pPredator:pPrey)
    //                );
    //    }
    //    BinomialDistribution startStateDist(startWeights);

    IntSampleStatistics startStateDist(
            PredPreyAgent::domainSize(),
            [](int agentId, int count) {
                double p = (PredPreyAgent(agentId).type() == PredPreyAgent::PREDATOR ? pPredator : pPrey);
                switch (count) {
                    case 0: return 1.0 - p;
                    case 1: return p;
                }
                return 0.0;
            });

    std::cout << "Initial state distribution = " << startStateDist << std::endl;


    AssimilationWindow <PredPreyAgent> window(
            windowSize,
            startStateDist,
            AgentStateObservation<PredPreyAgent>(
                    State<PredPreyAgent>(1,PredPreyAgent(0, 0, PredPreyAgent::PREY)),
                    1,
                    1.0
            )
    );

    std::cout << "Prior support is \n" << window.priorPMF.convexSupport << std::endl;
    std::cout << "Likelihood support is \n" << window.likelihoodPMF.convexSupport << std::endl;
    std::cout << "Posterior support is \n" << window.posterior.convexSupport << std::endl;

    MCMCSolver<PredPreyAgent> solver(window);
    solver.solve(5000,1000000);
    std::cout << "Analysis histograms:\n" << solver.solution << std::endl;
    std::cout << "Analysis = " << solver.solution.means() << std::endl;

    ExactSolver<CatMouseAgent> exact(window.posterior);
    std::cout << "Exact solution = " << exact.solution << std::endl;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void Experiments::CatMouseExpt() {
    AssimilationWindow<CatMouseAgent> window(
            2,
            BinomialDistribution({
                                         boost::math::binomial_distribution<double>(1.0, 0.9),
                                         boost::math::binomial_distribution<double>(1.0, 0.1),
                                         boost::math::binomial_distribution<double>(1.0, 0.1),
                                         boost::math::binomial_distribution<double>(1.0, 0.9)
                                 }),
            AgentStateObservation<CatMouseAgent>(
                    State<CatMouseAgent>(
                            1,
                            CatMouseAgent(CatMouseAgent::CAT, CatMouseAgent::LEFT)
                    ),
                    1,
                    1.0)
    );

    std::cout << "Prior support is \n" << window.priorPMF.convexSupport << std::endl;
    std::cout << "Likelihood support is \n" << window.likelihoodPMF.convexSupport << std::endl;
    std::cout << "Posterior support is \n" << window.posterior.convexSupport << std::endl;


    MCMCSolver<CatMouseAgent> solver(window);
    solver.solve(5000,1000000);
    std::cout << "Analysis = " << solver.solution.means() << std::endl;
//    window.doAnalysis(250000, 1000);
//    std::cout << window.analysis << std::endl;

    ExactSolver<CatMouseAgent> exact(window.posterior);
    std::cout << "Exact solution = " << exact.solution << std::endl;

}


void Experiments::CatMouseAssimilation() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    constexpr int windowSize = 2;
    constexpr double pMakeObservation = 0.25;//0.04;    // prob of making an observation of each gridsquare at each timestep
    constexpr double pObserveIfPresent = 1.0;
    constexpr int nSamplesPerWindow = 1000000; //250000;
    constexpr int nBurninSamples = 5000;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////

    BinomialDistribution startStateDist({
        boost::math::binomial_distribution<double>(1.0, 0.9),
        boost::math::binomial_distribution<double>(1.0, 0.1),
        boost::math::binomial_distribution<double>(1.0, 0.1),
        boost::math::binomial_distribution<double>(1.0, 0.9)
    });

    DataAssimilation<CatMouseAgent> assimilation(startStateDist, pMakeObservation, pObserveIfPresent);
    const AssimilationWindow<CatMouseAgent> &window = assimilation.addWindow(windowSize, nBurninSamples, nSamplesPerWindow);

    std::cout << "Prior support is \n" << window.priorPMF.convexSupport << std::endl;
    std::cout << "Likelihood support is \n" << window.likelihoodPMF.convexSupport << std::endl;
    std::cout << "Posterior support is \n" << window.posterior.convexSupport << std::endl;

//    MCMCSolver<CatMouseAgent> solver(window);
//    solver.solve(nBurninSamples, nSamplesPerWindow);
//    std::cout << "Analysis = " << solver.solution.means() << std::endl;

    std::cout << "Analysis histograms:\n" << assimilation.analysis << std::endl;

    std::cout << "Analysis means = " << assimilation.analysis.means() << std::endl;

    ExactSolver<CatMouseAgent> exact(window.posterior);
    std::cout << "Exact solution = " << exact.solution << std::endl;

}


void Experiments::BinomialAgentAssimilation() {
    BinomialAgent::GRIDSIZE = 3;
//    BinomialDistribution startState({
//        boost::math::binomial_distribution(1.0, 1.0),
//        boost::math::binomial_distribution(1.0, 0.1),
//        boost::math::binomial_distribution(1.0, 0.0)
//    });

    // gets stuck at kappaRow = -1
    BinomialDistribution startState({
        boost::math::binomial_distribution(1.0, 1.0),
        boost::math::binomial_distribution(1.0, 0.1),
        boost::math::binomial_distribution(1.0, 0.1)
    });


    AgentStateObservation<BinomialAgent> observation(State<BinomialAgent>(1, 0),1,0.9);

    AssimilationWindow<BinomialAgent> window(2, startState, observation);

    MCMCSolver<BinomialAgent> solver(window);
    solver.solve(5000, 1000000);
    std::cout << solver.solution.means() << std::endl;

    ExactSolver<BinomialAgent> exactSolver(window.posterior);
    std::cout << exactSolver.solution << std::endl;

}

void Experiments::RandomWalk() {
    using glp::X;
    glp::Problem myProb;

    myProb.addConstraint(1.0*X(1) + 1.0*X(2) + 1.0*X(3) <= 100.0);
    myProb.addConstraint(10.0*X(1) + 4.0*X(2) + 5.0*X(3) <= 600.0);
    myProb.addConstraint(2.0*X(1) + 2.0*X(2) + 6.0*X(3) <= 300.0);
    myProb.addConstraint(0.0 <= 1.0*X(1));
    myProb.addConstraint(0.0 <= 1.0*X(2));
    myProb.addConstraint(0.0 <= 1.0*X(3));
    myProb.setObjective(10.0*X(1) + 6.0*X(2) + 4.0*X(3));
    myProb.stdBasis();
    myProb.warmUp();
    std::cout << myProb;

    SimplexMCMC myMCMC(myProb, nullPMF);

    std::cout << myMCMC << std::endl;

    for(int sample=0; sample < 5; ++sample) {
        myMCMC.randomWalk();
        std::cout << myMCMC << std::endl;
        std::cout << "Sample is: " << myMCMC.X() << std::endl;
    }

//    mySimplex.pivot(3,3);
//    std::cout << mySimplex << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

double Experiments::informationGain(const std::vector<double> &realEndState, const IntSampleStatistics &prior, const IntSampleStatistics &analysis) {
    assert(realEndState.size() == prior.nDimensions());
    assert(prior.nDimensions() == analysis.nDimensions());
    std::vector<double> priorMeans = prior.means();
    std::vector<double> analysisMeans = analysis.means();
    for(int i=0; i<realEndState.size(); ++i) {
        std::cout << "Real occupancy = " << realEndState[i]
        << " prior p = " << prior.P(i, realEndState[i])
        << " posterior p = " << analysis.P(i, realEndState[i])
        << " post / prior p = " << analysis.P(i, realEndState[i])/prior.P(i, realEndState[i])
        << std::endl;
    }
    return (analysis.logP(realEndState) - prior.logP(realEndState))/log(2.0);
}


//void Experiments::GnuplotTest() {
//    Gnuplot gp; //(stdout);
//    auto state = ModelState<PredPreyAgent>::randomPoissonState([](const PredPreyAgent &agent) {
//        if(agent.type() == PredPreyAgent::PREDATOR) return 0.08;
//        return 0.16;
//    });
//
//    plotHeatMap(gp, state, state);
//}

//Gnuplot &Experiments::plotHeatMap(Gnuplot &gp, const PoissonState<PredPreyAgent> &aggregateState,
//                                  const ModelState<PredPreyAgent> &realState) {
//    typedef std::tuple<double,double,double,double,double> HeatRecord;
//    std::vector<std::vector<HeatRecord>> heatData;
//    std::vector<std::tuple<double,double,double>> pointData;
//
////    for(auto [agent, occupancy] : realState) {
////        pointData.emplace_back(agent.xPosition(), agent.yPosition(), agent.type()==PredPreyAgent::PREY?1:2);
////    }
//
//    for(int x=0; x<PredPreyAgent::GRIDSIZE; ++x) {
//        for(int y=0; y<PredPreyAgent::GRIDSIZE; ++y) {
//            int colour = 2*(realState[PredPreyAgent(x,y,PredPreyAgent::PREDATOR)]>0.0)
//                    + (realState[PredPreyAgent(x,y,PredPreyAgent::PREY)]>0.0);
//            if(colour != 0)
//                pointData.emplace_back(x, y, colour);
//        }
//    }
//
//
////    double predMaxOccupancy = 0.0;
////    double preyMaxOccupancy = 0.0;
////    for(auto [agent, occupancy] : aggregateState) {
////        if(agent.type() == PredPreyAgent::PREDATOR) {
////            if(occupancy > predMaxOccupancy) predMaxOccupancy = occupancy;
////        } else {
////            if(occupancy > preyMaxOccupancy) preyMaxOccupancy = occupancy;
////        }
////    }
//    double maxOccupancy = 0.0;
//    for(double occupancy : aggregateState.stateCounts) {
//        if(occupancy > maxOccupancy) maxOccupancy = occupancy;
//    }
////    std::cout << "Max occupancy = " << maxOccupancy << " at " << maxState << std::endl;
//
////    double predScale = 200.0/log(predMaxOccupancy + 1.0);
////    double preyScale = 200.0/log(preyMaxOccupancy + 1.0);
//    double predScale = 200.0/maxOccupancy; // predMaxOccupancy;
//    double preyScale = 200.0/maxOccupancy; // preyMaxOccupancy;
//    for(int x=0; x<PredPreyAgent::GRIDSIZE; ++x) {
//        std::vector<HeatRecord> &record = heatData.emplace_back();
//        for(int y=0; y<PredPreyAgent::GRIDSIZE; ++y) {
//            double nPrey = aggregateState.stateCounts[PredPreyAgent(x, y, PredPreyAgent::PREY)];
//            double nPred = aggregateState.stateCounts[PredPreyAgent(x, y, PredPreyAgent::PREDATOR)];
////            record.emplace_back(x, y, log(nPrey + 1.0) * preyScale, 0.0, log(nPred + 1.0) * predScale);
//            record.emplace_back(x,y,nPrey*preyScale,0.0,nPred*predScale);
//        }
//    }
//
//    gp << "set linetype 1 lc 'red'\n";
//    gp << "set linetype 2 lc 'blue'\n";
//    gp << "set linetype 3 lc 'magenta'\n";
//    gp << "plot [-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "] ";
//    gp << "'-' with rgbimage notitle, ";
//    gp << "'-' with points pointtype 5 pointsize 0.5 lc variable notitle\n";
//    gp.send2d(heatData);
//    gp.send1d(pointData);
//    return gp;
//}

//Gnuplot &Experiments::plotHeatMap(Gnuplot &gp, const BinomialDistribution &aggregateState,
//                                  const ModelState<PredPreyAgent> &realState) {
//    typedef std::tuple<double,double,double,double,double> HeatRecord;
//    std::vector<std::vector<HeatRecord>> heatData;
//    std::vector<std::tuple<double,double,double>> pointData;
//
//    for(int x=0; x<PredPreyAgent::GRIDSIZE; ++x) {
//        for(int y=0; y<PredPreyAgent::GRIDSIZE; ++y) {
//            int colour = 2*(realState[PredPreyAgent(x,y,PredPreyAgent::PREDATOR)]>0.0)
//                    + (realState[PredPreyAgent(x,y,PredPreyAgent::PREY)]>0.0);
//            if(colour != 0)
//                pointData.emplace_back(x, y, colour);
//        }
//    }
//
//    double maxOccupancy = 0.0;
//    for(double occupancy : aggregateState) {
//        if(occupancy > maxOccupancy) maxOccupancy = occupancy;
//    }
//    double predScale = 200.0/maxOccupancy; // predMaxOccupancy;
//    double preyScale = 200.0/maxOccupancy; // preyMaxOccupancy;
//    for(int x=0; x<PredPreyAgent::GRIDSIZE; ++x) {
//        std::vector<HeatRecord> &record = heatData.emplace_back();
//        for(int y=0; y<PredPreyAgent::GRIDSIZE; ++y) {
//            double nPrey = aggregateState.stateCounts[PredPreyAgent(x, y, PredPreyAgent::PREY)];
//            double nPred = aggregateState.stateCounts[PredPreyAgent(x, y, PredPreyAgent::PREDATOR)];
//            //            record.emplace_back(x, y, log(nPrey + 1.0) * preyScale, 0.0, log(nPred + 1.0) * predScale);
//            record.emplace_back(x,y,nPrey*preyScale,0.0,nPred*predScale);
//        }
//    }
//
//    gp << "set linetype 1 lc 'red'\n";
//    gp << "set linetype 2 lc 'blue'\n";
//    gp << "set linetype 3 lc 'magenta'\n";
//    gp << "plot [-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "] ";
//    gp << "'-' with rgbimage notitle, ";
//    gp << "'-' with points pointtype 5 pointsize 0.5 lc variable notitle\n";
//    gp.send2d(heatData);
//    gp.send1d(pointData);
//    return gp;
//}


Gnuplot &Experiments::plotAgents(Gnuplot &gp, const ModelState<PredPreyAgent> &state) {
    std::vector<std::tuple<double,double,double>> pointData;
    for(int x=0; x<PredPreyAgent::GRIDSIZE; ++x) {
        for(int y=0; y<PredPreyAgent::GRIDSIZE; ++y) {
            int colour = 2*(state[PredPreyAgent(x,y,PredPreyAgent::PREDATOR)]>0.0)
                         + (state[PredPreyAgent(x,y,PredPreyAgent::PREY)]>0.0);
            if(colour != 0)
                pointData.emplace_back(x, y, colour);
        }
    }

    gp << "set linetype 1 lc 'red'\n";
    gp << "set linetype 2 lc 'blue'\n";
    gp << "plot [-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "] ";
    gp << "'-' with points pointtype 5 pointsize 0.5 lc variable\n";
    gp.send1d(pointData);
    return gp;
}

