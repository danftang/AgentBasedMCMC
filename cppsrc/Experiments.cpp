//
// Created by daniel on 06/05/2021.
//

#include <vector>
#include "compose.h"
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
#include "BernoulliModelState.h"
#include "MCMCSampler.h"
#include "ModelStateSampleStatistics.h"
#include "RejectionSampler.h"

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

//    std::vector<boost::math::binomial_distribution<double>> startWeights(PredPreyAgent::domainSize());
//    for(int agentId=0; agentId < PredPreyAgent::domainSize(); ++agentId) {
//        startWeights[agentId] = boost::math::binomial_distribution<double>(
//                1,(PredPreyAgent(agentId).type() == PredPreyAgent::PREDATOR?pPredator:pPrey)
//                );
//    }

    BernoulliModelState<PredPreyAgent> startStateDist([pPredator,pPrey](PredPreyAgent agent) {
        return agent.type() == PredPreyAgent::PREDATOR?pPredator:pPrey;
    });

    AssimilationWindow<PredPreyAgent> window(windowSize, startStateDist, pMakeObservation, pObserveIfPresent);

    MCMCSampler sampler(window.posterior);
//    solver.solve(nBurnInSamples, nSamplesPerWindow);


    return std::vector<double>(); //assimilation.calculateInformationGain();
}


void Experiments::PredPreyAssimilation() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    PredPreyAgent::GRIDSIZE = 8;
    constexpr int windowSize = 2;
    constexpr int nWindows = 1;
    constexpr double pPredator = 0.08;//0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
    constexpr double pMakeObservation = 0.2;//0.04;    // prob of making an observation of each gridsquare at each timestep
    constexpr double pObserveIfPresent = 0.9;
    constexpr int nSamplesPerWindow = 1500000; //250000;
    constexpr int nBurninSamples = 1000;
    constexpr int nPriorSamples = 100000;
//    constexpr int plotTimestep = nTimesteps-1;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////

    BernoulliModelState<PredPreyAgent> startStateDist([pPredator,pPrey](PredPreyAgent agent) {
        return agent.type() == PredPreyAgent::PREDATOR?pPredator:pPrey;
    });
    std::cout << "Initial state distribution = " << startStateDist << std::endl;
    Trajectory<PredPreyAgent> realTrajectory(nWindows*windowSize, startStateDist.sampler());

    const Distribution<ModelState<PredPreyAgent>> *analysis = &startStateDist;
    ModelStateSampleStatistics<PredPreyAgent> sampleStats;
    for(int w=0; w<nWindows; ++w) {
        AssimilationWindow<PredPreyAgent> window(
                *analysis,
                realTrajectory.slice(w*windowSize, windowSize),
                pMakeObservation,
                pObserveIfPresent
        );
        MCMCSampler sampler(window.posterior, window.priorSampler());
        for(int s=0; s<nBurninSamples; ++s) sampler.nextSample();
        sampleStats.clear();
        sampleStats.sampleFromEndState(sampler, nSamplesPerWindow);
        analysis = &sampleStats;

        std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
        std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
        std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;
        std::cout << "Analysis means = " << sampleStats.means() << std::endl;

        ABMPlotter<PredPreyAgent> gp;
        gp.plot(window.realTrajectory.endState(), sampleStats.means());
//        ModelStateSampleStatistics priorEnd(window.priorSampler, 10000);
//        std::cout << "Window information gain = "
//        << informationGain(window.realTrajectory.endState(), priorEnd, sampleStats) << std::endl;
    }

    std::cout << "Calculating prior end state" << std::endl;
    ModelStateSampleStatistics<PredPreyAgent> priorEnd;
    priorEnd.sampleFromEndState(
            Trajectory<PredPreyAgent>::priorSampler(nWindows * windowSize, startStateDist.sampler()),
            nPriorSamples
    );

    std::cout << "Prior end means = " << priorEnd.means() << std::endl;
    std::cout << "Total information gain = " << informationGain(realTrajectory.endState(), priorEnd, sampleStats) << std::endl;
}


void Experiments::PredPreySingleObservation() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    PredPreyAgent::GRIDSIZE = 3;
    constexpr int nTimesteps = 2;
    constexpr double pPredator = 0.1;//0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
    constexpr int nSamples = 1500000; //250000;
    constexpr int nBurninSamples = 10000;
    constexpr int nRejectionSamples = 250000;
    //    constexpr int plotTimestep = nTimesteps-1;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////

    BernoulliModelState<PredPreyAgent> startStateDist([pPredator,pPrey](PredPreyAgent agent) {
        return agent.type() == PredPreyAgent::PREDATOR?pPredator:pPrey;
    });

    std::cout << "Initial state distribution = " << startStateDist << std::endl;

    AssimilationWindow <PredPreyAgent> window(
            nTimesteps,
            startStateDist,
            AgentStateObservation<PredPreyAgent>(
                    State<PredPreyAgent>(1,PredPreyAgent(1, 1, PredPreyAgent::PREY)),
                    1,
                    1.0
            )
    );

    std::cout << "Prior support is \n" << window.priorPMF.convexSupport << std::endl;
    std::cout << "Likelihood support is \n" << window.likelihoodPMF.convexSupport << std::endl;
    std::cout << "Posterior support is \n" << window.posterior.convexSupport << std::endl;
    std::cout << "Real trajectory is " << window.realTrajectory << std::endl;

    MCMCSampler sampler(window.posterior, window.priorSampler());
    for(int s=0; s<nBurninSamples; ++s) sampler.nextSample();
    ModelStateSampleStatistics<PredPreyAgent> sampleStats(sampler, nSamples);

    std::cout << "Analysis histograms:\n" << sampleStats << std::endl;
    std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
    std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
    std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;
    std::cout << "Analysis means = " << sampleStats.means() << std::endl;

    RejectionSampler rejectionSampler(window.priorSampler, window.likelihoodPMF);
    ModelStateSampleStatistics<PredPreyAgent> rejectionSampleStats(rejectionSampler, nRejectionSamples);
    std::cout << "Rejection means = " << rejectionSampleStats.means() << std::endl;

    ABMPlotter<PredPreyAgent> gp;
    gp.plot(window.realTrajectory.endState(), sampleStats.means());

    ABMPlotter<PredPreyAgent> gp2;
    gp2.plot(window.realTrajectory.endState(), rejectionSampleStats.means());

}


void Experiments::CatMouseAssimilation() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    constexpr int nWindows = 2;
    constexpr int windowSize = 2;
    constexpr double pMakeObservation = 0.3;//0.04;    // prob of making an observation of each gridsquare at each timestep
    constexpr double pObserveIfPresent = 1.0;
    constexpr int nSamplesPerWindow = 1500000; //250000;
    constexpr int nBurninSamples = 5000;
    constexpr int nPriorSamples = 100000;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////

    BernoulliModelState<CatMouseAgent> startStateDist({0.9, 0.1, 0.1, 0.9});

    Trajectory<CatMouseAgent> realTrajectory(nWindows*windowSize, startStateDist.sampler());
    const Distribution<ModelState<CatMouseAgent>> *analysis = &startStateDist;
    ModelStateSampleStatistics<CatMouseAgent> sampleStats;
    for(int w=0; w<nWindows; ++w) {
//        AssimilationWindow<CatMouseAgent> window(windowSize, analysis, pMakeObservation, pObserveIfPresent);
        AssimilationWindow<CatMouseAgent> window(
                *analysis,
                realTrajectory.slice(w*windowSize, windowSize),
                pMakeObservation,
                pObserveIfPresent
        );
        std::cout << "Prior support is \n" << window.priorPMF.convexSupport << std::endl;
        std::cout << "Likelihood support is \n" << window.likelihoodPMF.convexSupport << std::endl;
        std::cout << "Posterior support is \n" << window.posterior.convexSupport << std::endl;

        MCMCSampler sampler(window.posterior, window.priorSampler());
        for(int s=0; s<nBurninSamples; ++s) sampler.nextSample();
        sampleStats.clear();
        sampleStats.sampleFromEndState(sampler, nSamplesPerWindow);
        analysis = &sampleStats;

        std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
        std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
        std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;

        std::cout << "Analysis histograms =\n" << sampleStats << std::endl;
        std::cout << "Analysis Means = " << sampleStats.means() << std::endl;

        ExactSolver<CatMouseAgent> exact(window.posterior);
        std::cout << "Exact exactEndState = " << exact.exactEndState << std::endl;
    }

    ModelStateSampleStatistics<CatMouseAgent> priorEnd;
    priorEnd.sampleFromEndState(
            Trajectory<CatMouseAgent>::priorSampler(nWindows * windowSize, startStateDist.sampler()),
            nPriorSamples
    );
    std::cout << "Prior end means = " << priorEnd.means() << std::endl;
    std::cout << "Real end state = " << realTrajectory.endState() << std::endl;
    std::cout << "Total information gain = " << informationGain(realTrajectory.endState(), priorEnd, sampleStats) << std::endl;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void Experiments::CatMouseSingleObservation() {
    // TODO: analysis probs slightly out at 3 timesteps (with infeasibleExpectationFraction = 0.5. Closer when 1.0)
    constexpr int nTimesteps = 3;
    constexpr int nBurninSamples = 10000;
    constexpr int nSamples = 1000000;

    AssimilationWindow<CatMouseAgent> window(
            nTimesteps,
            BernoulliModelState<CatMouseAgent>({0.9, 0.1, 0.1, 0.9}),
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

    MCMCSampler sampler(window.posterior, window.priorSampler());
    for(int s=0; s<nBurninSamples; ++s) sampler.nextSample();
    ModelStateSampleStatistics<CatMouseAgent> sampleStats(sampler, nSamples);

    std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
    std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
    std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;
    std::cout << "Analysis means = " << sampleStats.means() << std::endl;

    ExactSolver<CatMouseAgent> exact(window.posterior);
    std::cout << "Exact means = " << exact.exactEndState << std::endl;

}


void Experiments::BinomialAgentAssimilation() {
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

    MCMCSampler<Trajectory<BinomialAgent>> sampler(window.posterior, window.priorSampler());
    std::cout << "simplex = \n" << sampler.simplex << std::endl;
    std::cout << "Starting burn-in" << std::endl;
//    for(int s=0; s<nBurninSamples; ++s) sampler();
    std::cout << "Done burn-in" << std::endl;
    ModelStateSampleStatistics<BinomialAgent> sampleStats(sampler, nSamplesPerWindow);

    std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
    std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
    std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;
    std::cout << "Analysis means = " << sampleStats.means() << std::endl;

//    std::cout << "Starting exact solve" << std::endl;
    ExactSolver<BinomialAgent> exactSolver(window.posterior);
    std::cout << "Exact means = " << exactSolver.exactEndState << std::endl;

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
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


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

