//
// Created by daniel on 06/05/2021.
//

#include <vector>
#include <thread>
#include <future>
#include "compose.h"
#include "Experiments.h"
#include "agents/CatMouseAgent.h"
#include "PredPreyProblem.h"
#include "SimplexMCMC.h"
#include "agents/PredPreyAgent.h"
#include "PoissonState.h"
#include "DataAssimilation.h"
#include "debug.h"
#include "BinomialDistribution.h"
#include "SampleStatistics.h"
#include "TrajectorySampler.h"
#include "TrajectoryLikelihoodPMF.h"
#include "ExactSolver.h"
#include "agents/BinomialAgent.h"
#include "MCMCSolver.h"
#include "Plotter.h"
#include "StlStream.h"
#include "BernoulliModelState.h"
#include "MCMCSampler.h"
#include "ModelStateSampleStatistics.h"
#include "RejectionSampler.h"
#include "diagnostics/MeanAndVariance.h"
#include "diagnostics/ValarrayLogger.h"
#include "diagnostics/ScaleReduction.h"
#include "diagnostics/Autocorrelation.h"
#include "GnuplotExtensions.h"
#include "diagnostics/Dataflow.h"
#include "diagnostics/MultiChainStats.h"


void Experiments::animatedPredPreyDemo() {
    constexpr int GRIDSIZE=64;
    int nTimesteps = 128;
    int framerate = 15;

    std::cout << PredPreyAgent<GRIDSIZE>::pPredBirthGivenPrey
              << " " << PredPreyAgent<GRIDSIZE>::pPredDie
              << " " << PredPreyAgent<GRIDSIZE>::pPreyEatenGivenPred
              << " " << PredPreyAgent<GRIDSIZE>::pPreyDie
              << " " << PredPreyAgent<GRIDSIZE>::pPreyBirth
              << " " << std::endl;

    auto startStatePMF = BernoulliModelState<PredPreyAgent<GRIDSIZE>>([](PredPreyAgent<GRIDSIZE> agent) {
        return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?0.05:0.05;
    });

    Random::seedFromTimeNow();
    Trajectory<PredPreyAgent<GRIDSIZE>> traj(nTimesteps, startStatePMF.sampler());
    Plotter().animate(traj, framerate);
}


void Experiments::DataflowDemo() {
    using namespace dataflow;
    int nBurnIn = 100;
    int nSamples = 200;
    auto sampler = [n = 0]() mutable { return n++; };
    std::function<bool(int)> cint = [](int x) { std::cout << x << std::endl; return true; };
    auto cany = [](auto x) { std::cout << x << std::endl; return true; };

    auto trajectoryToEnergy = [](int i) { return -i*i; };
    auto synopsis = [](int i) { return std::vector<double>{ 1.0*i, 2.0*i }; };
    MeanAndVariance  meanVariances;
    std::vector<std::vector<double>> measureLog;

//    sampler >>= Drop{nBurnIn} >>= Split {
//            Thin(10) >>= Map { trajectoryToEnergy } >>= plot1DAfter(nSamples/10, "Energy"),
//            Take(nSamples) >>= Map { synopsis } >>= Split{
//                    meanVariances.consumer(),
//                    save(measureLog)
//            }
//    };

    std::valarray<int> log1(5);
    std::valarray<int> log2(4);
    std::valarray<int> log3(1);


    sampler >>= SwitchOnClose {
        save(log1),
        save(log2),
        save(log3)
    };

    std::cout << log1 << std::endl;
    std::cout << log2 << std::endl;
    std::cout << log3 << std::endl;

//    std::cout << meanVariances.mean() << std::endl;
//    std::cout << meanVariances.sampleVariance() << std::endl;
//    std::cout << measureLog << std::endl;


}

template<int GRIDSIZE>
MultiChainStats Experiments::PredPreyConvergenceThread(const ConvexPMF<Trajectory<PredPreyAgent<GRIDSIZE>>> &posterior, Trajectory<PredPreyAgent<GRIDSIZE>> startState) {
    using namespace dataflow;
    constexpr int nSamples = 50000; // must be an even number
    assert((nSamples&1) == 0);
    const double maxLagProportion = 0.4;
    const int nLags = 80;
    constexpr int nBurnIn = nSamples*0.5;
    auto trajectoryToEnergy = [](const Trajectory<PredPreyAgent<GRIDSIZE>> &trajectory) { return -trajectory.logProb(); };
    auto trajectoryToEndState = [](const Trajectory<PredPreyAgent<GRIDSIZE>> &trajectory) { return trajectory.endState(); };
    MCMCSampler sampler(posterior, startState);

    std::valarray<std::valarray<double>> firstSynopsisSamples(nSamples/2);
    std::valarray<std::valarray<double>> lastSynopsisSamples(nSamples/2);
    ModelState<PredPreyAgent<GRIDSIZE>> firstMeanEndState;
    ModelState<PredPreyAgent<GRIDSIZE>> lastMeanEndState;
    std::valarray<Trajectory<PredPreyAgent<GRIDSIZE>>> firstNextSample(Trajectory<PredPreyAgent<GRIDSIZE>>(0),1);
    std::valarray<Trajectory<PredPreyAgent<GRIDSIZE>>> lastNextSample(Trajectory<PredPreyAgent<GRIDSIZE>>(0),1);


//    sampler >>= Split {
//            Thin(10) >>= Map { trajectoryToEnergy } >>= plot1DAfter(nSamples/10, "Energy"),
//            Drop(nBurnIn) >>= Take((nSamples/2)*2) >>= Map { Experiments::Synopsis } >>= SwitchAfter {nSamples / 2,
//                    save(firstSynopsisSamples),
//                    save(lastSynopsisSamples)
//            }
//    };

    sampler >>= Drop(nBurnIn)
            >>= SwitchOnClose {
                    Split {
                        Map{ Experiments::Synopsis<GRIDSIZE> } >>= save(firstSynopsisSamples),
                        Take(nSamples/2) >>= Map{ trajectoryToEndState } >>= Sum(firstMeanEndState)
                    },
                    save(firstNextSample),
                    Split{
                        Map{ Experiments::Synopsis<GRIDSIZE> } >>= save(lastSynopsisSamples),
                        Take(nSamples/2) >>= Map{ trajectoryToEndState } >>= Sum(firstMeanEndState)
                    },
                    save(lastNextSample)
                };


    std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
    std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
    std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;

//    std::cout << firstSynopsisSamples << std::endl;
//    std::cout << lastSynopsisSamples << std::endl;

    firstMeanEndState /= nSamples/2;
    lastMeanEndState /= nSamples/2;

    MultiChainStats stats;
    stats.reserve(2);
    stats.emplace_back(
            std::move(firstSynopsisSamples),
            nLags,
            maxLagProportion,
            std::move(firstMeanEndState),
            std::move(firstNextSample[0]),
            sampler.simplex.feasibleStatistics,
            sampler.simplex.infeasibleStatistics);
    stats.emplace_back(
            std::move(lastSynopsisSamples),
            nLags,
            maxLagProportion,
            std::move(lastMeanEndState),
            std::move(lastNextSample[0]),
            sampler.simplex.feasibleStatistics,
            sampler.simplex.infeasibleStatistics);
//    std::cout << stats << std::endl;
    return std::move(stats);
}


void Experiments::PredPreyConvergence() {
    constexpr int GRIDSIZE = 8;
    constexpr int windowSize = 2;
    constexpr double pPredator = 0.05;//0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
    constexpr double pMakeObservation = 0.1;//0.04;    // prob of making an observation of each gridsquare at each timestep
    constexpr double pObserveIfPresent = 0.9;
    constexpr int nThreads = 4;

    std::stringstream strstr;
    boost::archive::binary_oarchive boostout(strstr);

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////

    BernoulliModelState<PredPreyAgent<GRIDSIZE>> startStateDist([pPredator,pPrey](PredPreyAgent<GRIDSIZE> agent) {
        return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
    });
    std::cout << "Initial state distribution = " << startStateDist << std::endl;
    Trajectory<PredPreyAgent<GRIDSIZE>> realTrajectory(windowSize, startStateDist.sampler());

    ModelStateSampleStatistics<PredPreyAgent<GRIDSIZE>> sampleStats;
    AssimilationWindow<PredPreyAgent<GRIDSIZE>> window(startStateDist,realTrajectory, pMakeObservation, pObserveIfPresent);

    std::future<MultiChainStats> futureResults[nThreads];
    for(int thread = 0; thread < nThreads; ++thread) {
        futureResults[thread] = std::async(&PredPreyConvergenceThread<GRIDSIZE>, window.posterior, window.priorSampler());
    }

//    std::vector<MeanAndVariance> meanvariances;
    MultiChainStats multiChainStats;
    multiChainStats.reserve(2*nThreads);
    for(int thread=0; thread<nThreads; ++thread) {
        futureResults[thread].wait();
        multiChainStats += futureResults[thread].get();
    }

    boostout << multiChainStats;

    MultiChainStats statsReloaded;
    boost::archive::binary_iarchive boostin(strstr);
    boostin >> statsReloaded;

    std::cout << statsReloaded;

    Plotter gp;
    std::valarray<std::valarray<double>> autocorrelation = multiChainStats.autocorrelation();
    gp.heatmap(autocorrelation, 0.5, autocorrelation[0].size()-0.5, 0.5*multiChainStats.front().varioStride, (autocorrelation.size()+0.5)*multiChainStats.front().varioStride);

    Plotter gp2;
    gp2 << "plot '-' using 0:1 with lines, '-' using 0:2 with lines, 0 with lines\n";
    gp2.send1d(autocorrelation);
    gp2.send1d(autocorrelation);

    std::valarray<double> neff = multiChainStats.effectiveSamples();
    std::valarray<double> ineff = (multiChainStats.nSamples()*1.0)/neff;

    std::cout << multiChainStats << std::endl;
    std::cout << "Potential scale reduction: " << multiChainStats.potentialScaleReduction() << std::endl;
    std::cout << "Actual number of samples per chain: " << multiChainStats.nSamples() << std::endl;
    std::cout << "Effective number of samples: " << neff << std::endl;
    std::cout << "Sample inefficiency factor: " << ineff << std::endl;

    Plotter gp3;
    gp3.plot(realTrajectory.endState(), multiChainStats.meanEndState());

//    std::valarray<double> gelman = gelmanScaleReduction(meanvariances);
//    std::cout << std::endl << "Gelman scale reduction factors: " << gelman << std::endl;
}


void Experiments::PredPreyAssimilation() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    // 8 x 8 x 8 45% infeasible, 400 iterations per transtition, 20ms per transition
    constexpr int GRIDSIZE = 16;
    constexpr int windowSize = 8;
    constexpr double pPredator = 0.05;//0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 0.05;    // Poisson prob of prey in each gridsquare at t=0
    constexpr double pMakeObservation = 0.02;//0.04;    // prob of making an observation of each gridsquare at each timestep
    constexpr double pObserveIfPresent = 0.9;
    constexpr int nSamplesPerWindow = 50000; //250000;
    constexpr int nBurninSamples = 5000;
    constexpr int nPriorSamples = 100000;
//    constexpr int plotTimestep = nTimesteps-1;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////

    PoissonModelState<PredPreyAgent<GRIDSIZE>> startStateDist([](PredPreyAgent<GRIDSIZE> agent) {
        return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
    });

//    BernoulliModelState<PredPreyAgent<GRIDSIZE>> startStateDist([](PredPreyAgent<GRIDSIZE> agent) {
//        return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
//    });

    std::cout << "Initial state distribution = " << startStateDist << std::endl;
    Trajectory<PredPreyAgent<GRIDSIZE>> realTrajectory(windowSize, startStateDist.sampler());

//    const Distribution<ModelState<PredPreyAgent<GRIDSIZE>>> *analysis = &startStateDist;
    ModelStateSampleStatistics<PredPreyAgent<GRIDSIZE>> sampleStats;
    AssimilationWindow<PredPreyAgent<GRIDSIZE>> window(
            startStateDist,
            realTrajectory,
            pMakeObservation,
            pObserveIfPresent);
    MCMCSampler sampler(window.posterior, window.priorSampler(), Trajectory<PredPreyAgent<GRIDSIZE>>::marginalLogProbsByEvent(windowSize));
    for(int s=0; s<nBurninSamples; ++s) sampler.nextSample();
    sampleStats.clear();
    auto startTime = std::chrono::high_resolution_clock::now();
    sampleStats.sampleFromEndState(sampler, nSamplesPerWindow);
    auto endTime = std::chrono::high_resolution_clock::now();
    auto execTime = endTime - startTime;
    std::cout << "Finished sampling in " << execTime << std::endl;

    std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
    std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
    std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/(sampler.simplex.feasibleStatistics.nSamples + sampler.simplex.infeasibleStatistics.nSamples) << "%" << std::endl;
    std::cout << "Iterations per feasible transition = " << (sampler.simplex.feasibleStatistics.nSamples + sampler.simplex.infeasibleStatistics.nSamples) * 1.0/sampler.simplex.feasibleStatistics.nNonDegenerate << std::endl;
    std::cout << "Milliseconds per feasible transition = " << execTime.count() * 1e-6/sampler.simplex.feasibleStatistics.nNonDegenerate << std::endl;
    std::cout << "Analysis means = " << sampleStats.means() << std::endl;


    Plotter gp;
    gp.plot(window.realTrajectory.endState(), sampleStats.means());
//        ModelStateSampleStatistics priorEnd(window.priorSampler, 10000);
//        std::cout << "Window information gain = "
//        << informationGain(window.realTrajectory.endState(), priorEnd, sampleStats) << std::endl;
//    std::cout << "Calculating prior end state" << std::endl;
//    ModelStateSampleStatistics<PredPreyAgent<GRIDSIZE>> priorEnd;
//    priorEnd.sampleFromEndState(
//            Trajectory<PredPreyAgent<GRIDSIZE>>::priorSampler(nWindows * windowSize, startStateDist.sampler()),
//            nPriorSamples
//    );
//
//    std::cout << "Prior end means = " << priorEnd.means() << std::endl;
//    std::cout << "Total information gain = " << informationGain(realTrajectory.endState(), priorEnd, sampleStats) << std::endl;
}


void Experiments::PredPreySingleObservation() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    constexpr int GRIDSIZE = 3;
    constexpr int nTimesteps = 2;
    constexpr double pPredator = 0.1;//0.08;          // Poisson prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Poisson prob of prey in each gridsquare at t=0
    constexpr int nSamples = 1500000; //250000;
    constexpr int nBurninSamples = 10000;
    constexpr int nRejectionSamples = 500000;
    //    constexpr int plotTimestep = nTimesteps-1;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////

    BernoulliModelState<PredPreyAgent<GRIDSIZE>> startStateDist([pPredator,pPrey](PredPreyAgent<GRIDSIZE> agent) {
        return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
    });

    std::cout << "Initial state distribution = " << startStateDist << std::endl;

    AssimilationWindow <PredPreyAgent<GRIDSIZE>> window(
            nTimesteps,
            startStateDist,
            AgentStateObservation<PredPreyAgent<GRIDSIZE>>(
                    State<PredPreyAgent<GRIDSIZE>>(1,PredPreyAgent<GRIDSIZE>(1, 1, PredPreyAgent<GRIDSIZE>::PREY)),
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
    ModelStateSampleStatistics<PredPreyAgent<GRIDSIZE>> sampleStats(sampler, nSamples);

    std::cout << "Analysis histograms:\n" << sampleStats << std::endl;
    std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
    std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
    std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;
    std::cout << "Analysis means = " << sampleStats.means() << std::endl;

    RejectionSampler rejectionSampler(window.priorSampler, window.likelihoodPMF);
    ModelStateSampleStatistics<PredPreyAgent<GRIDSIZE>> rejectionSampleStats(rejectionSampler, nRejectionSamples);
    std::cout << "Rejection means = " << rejectionSampleStats.means() << std::endl;

    Plotter gp;
    gp.plot(window.realTrajectory.endState(), sampleStats.means());

    Plotter gp2;
    gp2.plot(window.realTrajectory.endState(), rejectionSampleStats.means());

}


void Experiments::CatMouseAssimilation() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    constexpr int nWindows = 2;
    constexpr int windowSize = 2;
    constexpr double pMakeObservation = 0.3;//0.04;    // prob of making an observation of each gridsquare at each timestep
    constexpr double pObserveIfPresent = 1.0;
    constexpr int nSamplesPerWindow = 1500000; //250000;
    constexpr int nBurninSamples = 10000;
//    constexpr int nPriorSamples = 100000;

//    with all agent states at T=0 observed 1, 0 ,0 ,1
//    Analysis Means = {0.609056, 0.390944, 0.280905, 0.719095}
//    Exact exactEndState = {0.625, 0.375, 0.25, 0.75}
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

//    ModelStateSampleStatistics<CatMouseAgent> priorEnd;
//    priorEnd.sampleFromEndState(
//            Trajectory<CatMouseAgent>::priorSampler(nWindows * windowSize, startStateDist.sampler()),
//            nPriorSamples
//    );
//    std::cout << "Prior end means = " << priorEnd.means() << std::endl;
//    std::cout << "Real end state = " << realTrajectory.endState() << std::endl;
//    std::cout << "Total information gain = " << informationGain(realTrajectory.endState(), priorEnd, sampleStats) << std::endl;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void Experiments::CatMouseMultiObservation() {
    constexpr int nTimesteps = 2;
//    constexpr int nBurninSamples = 1000;
    constexpr int nSamples = 100000;

    Random::gen.seed(530673);

    BernoulliModelState<CatMouseAgent>      startState({0.9, 0.1, 0.1, 0.9});
    ConvexPMF<Trajectory<CatMouseAgent>>    prior(nTimesteps, startState.PMF());
    ConvexPMF<Trajectory<CatMouseAgent>> likelihood(nTimesteps, {
            AgentStateObservation<CatMouseAgent>(
                    State<CatMouseAgent>(0,CatMouseAgent(CatMouseAgent::CAT, CatMouseAgent::LEFT)),
                    1,
                    1.0
            ),
            AgentStateObservation<CatMouseAgent>(
                    State<CatMouseAgent>(0,CatMouseAgent(CatMouseAgent::CAT, CatMouseAgent::RIGHT)),
                    0,
                    1.0
            ),
            AgentStateObservation<CatMouseAgent>(
                    State<CatMouseAgent>(0,CatMouseAgent(CatMouseAgent::MOUSE, CatMouseAgent::LEFT)),
                    0,
                    1.0
            ),
            AgentStateObservation<CatMouseAgent>(
                    State<CatMouseAgent>(0,CatMouseAgent(CatMouseAgent::MOUSE, CatMouseAgent::RIGHT)),
                    1,
                    1.0
            )
    });
    ConvexPMF<Trajectory<CatMouseAgent>> posterior(prior*likelihood);

    std::cout << "Prior support is \n" << prior.convexSupport << std::endl;
    std::cout << "Likelihood support is \n" << likelihood.convexSupport << std::endl;
    std::cout << "Posterior support is \n" << posterior.convexSupport << std::endl;

    MCMCSampler sampler(posterior, Trajectory<CatMouseAgent>(
//            {0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1} // cat stay,stay
    {0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0} // cat move,stay
    ));
    std::cout << "Simplex geometry = " << sampler.simplex.nBasic() << " x " << sampler.simplex.nNonBasic() << std::endl;
    //    for(int s=0; s<nBurninSamples; ++s) sampler.nextSample();

    std::multiset<Trajectory<CatMouseAgent>> trajHistogram;
    for(int s=0; s<nSamples; ++s) {
        trajHistogram.insert(sampler.nextSample());
    }

//    ModelStateSampleStatistics<CatMouseAgent> sampleStats(sampler, nSamples);

    std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
    std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
    std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;

    //    std::cout << "Analysis means = " << sampleStats.means() << std::endl;
    double mse = 0.0;
    double exact[] = {0.5625, 0.1875, 0.1875, 0.0625};
    int i=0;
    for(auto Tp = trajHistogram.begin(); Tp != trajHistogram.end(); Tp = trajHistogram.upper_bound(*Tp)) {
        double pState = trajHistogram.count(*Tp)*1.0/nSamples;
        mse += pow(pState - exact[i],2.0);
        std::cout << pState << " " << *Tp << std::endl;
        ++i;
    }

    std::cout << "Exact values (calculated by hand) should be: 0.5625, 0.1875, 0.1875, 0.0625" << std::endl;
    std::cout << "Mean square error = " << mse << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void Experiments::CatMouseSingleObservation() {
    // TODO: analysis probs slightly out at 3 timesteps (with infeasibleExpectationFraction = 0.5. Closer when 1.0)
    constexpr int nTimesteps = 3;
    constexpr int nBurninSamples = 2000;
    constexpr int nSamples = 20000;

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

    constexpr int NCHAINS = 4;
    std::vector<std::vector<double>>                        energies(NCHAINS);
    std::vector<MeanAndVariance>                            meanVariances(NCHAINS,MeanAndVariance());
    std::vector<ValarrayLogger<std::valarray<double>>>      loggers(NCHAINS, ValarrayLogger<std::valarray<double>>(nSamples));
    std::vector<MCMCSampler<Trajectory<CatMouseAgent>>>     samplers;
    samplers.reserve(NCHAINS);

    for(int chain=0; chain<NCHAINS; ++chain) {
        std::cout << "Starting thread " << chain << std::endl;
        auto &sampler = samplers.emplace_back(window.posterior, window.priorSampler());
        sampler.burnIn(nBurninSamples);
        sampler.sampleInThread(nSamples,[&vaLogger = loggers[chain], &mvLogger=meanVariances[chain], &E = energies[chain]](const Trajectory<CatMouseAgent> &sample) {
            mvLogger(sample);
            vaLogger(std::valarray<double>(sample.data(), sample.size()));
            E.push_back(-sample.logProb());
        });
//        sampler.sample(nSamples,[&vaLogger = loggers[chain], &mvLogger=meanVariances[chain]](const Trajectory<CatMouseAgent> &sample) {
//            mvLogger(sample);
//            vaLogger(std::valarray<double>(sample.data(), sample.size()));
//        });
    }


//    MCMCSampler sampler(window.posterior, window.priorSampler());
//    for(int s=0; s<nBurninSamples; ++s) sampler.nextSample();
//    ModelStateSampleStatistics<CatMouseAgent> sampleStats(sampler, nSamples);


    for(int thread=0; thread<NCHAINS; ++thread) {
        samplers[thread].join();
        std::cout << "Statistics for chain number " << thread << std::endl;
        std::cout << "Feasible stats =\n" << samplers[thread].simplex.feasibleStatistics << std::endl;
        std::cout << "Infeasible stats =\n" << samplers[thread].simplex.infeasibleStatistics << std::endl;
        std::cout << "Infeasible proportion = "
                  << samplers[thread].simplex.infeasibleStatistics.nSamples * 100.0 / samplers[thread].simplex.feasibleStatistics.nSamples
                  << "%" << std::endl;
        std::valarray<std::valarray<double>> ac = geyerAutocorrelation(loggers[thread].samples, 40, 0.2);
        std::cout << "Geyer autocorrelation: " << ac << std::endl;
        Gnuplot gp;
        gp << ac;
        Gnuplot gp2;
        gp2 << "plot '-' with lines\n";
        gp2.send1d(energies[thread]);
    }
    std::valarray<double> gelman = gelmanScaleReduction(meanVariances);
    std::cout << "Gelman scale reduction: " << gelman << std::endl;
    std::vector<double> gelmanV(std::begin(gelman),std::end(gelman));
    Gnuplot gp;
    gp << "plot '-' with lines\n";
    gp.send1d(gelmanV);
//    std::cout << "Analysis means = " << sampleStats.means() << std::endl;

//    ExactSolver<CatMouseAgent> exact(window.posterior);
//    std::cout << "Exact means = " << exact.exactEndState << std::endl;

}


void Experiments::CatMousePrior() {
    ////////////////////////////////////////// SETUP PARAMETERS ////////////////////////////////////////
    constexpr int nTimesteps = 2;
    constexpr int nSamples = 1500000; //250000;
    constexpr int nBurninSamples = 5000;
    constexpr int nPriorSamples = 100000;

    ////////////////////////////////////////// SETUP PROBLEM ////////////////////////////////////////

    BernoulliModelState<CatMouseAgent> startStateDist({0.9, 0.1, 0.1, 0.9});
    ConvexPMF<Trajectory<CatMouseAgent>> priorPMF(nTimesteps, startStateDist.PMF());
    std::cout << "Prior support = \n" << priorPMF.convexSupport << std::endl;

    MCMCSampler sampler(priorPMF, Trajectory<CatMouseAgent>(nTimesteps, startStateDist.sampler()));
    for(int s=0; s<nBurninSamples; ++s) sampler.nextSample();
    ModelStateSampleStatistics<CatMouseAgent> sampleStats(sampler, nSamples);

    std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
    std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
    std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;

    std::cout << "Analysis histograms =\n" << sampleStats << std::endl;
    std::cout << "Analysis Means = " << sampleStats.means() << std::endl;

    ExactSolver<CatMouseAgent> exact(priorPMF);
    std::cout << "Exact exactEndState = " << exact.exactEndState << std::endl;
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
    for(int s=0; s<nBurninSamples; ++s)  sampler();
    std::cout << "Done burn-in" << std::endl;
    ModelStateSampleStatistics<BinomialAgent> sampleStats(sampler, nSamplesPerWindow);

    std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
    std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
    std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;
    std::vector<double> analysis = sampleStats.means();
    std::cout << "Analysis means = " << analysis << std::endl;

//    std::cout << "Starting exact solve" << std::endl;
    ExactSolver<BinomialAgent> exactSolver(window.posterior);
    std::cout << "Exact means = " << exactSolver.exactEndState << std::endl;
    double mse = 0.0;
    for(int i=0; i<analysis.size(); ++i) {
        mse += pow(analysis[i] - exactSolver.exactEndState[i],2.0);
    }
    mse /= analysis.size();
    std::cout << "Mean square error = " << mse << std::endl;
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


// Test whether all vertices of the polyhedron described by the
// Fermionic + continuity constraints are all integer
void Experiments::FermionicIntegrality() {
    constexpr int GRIDSIZE = 3;
    int nTimesteps = 2;
    ConvexPolyhedron fermi =
            ABMConstraints<PredPreyAgent<GRIDSIZE>>::actFermionicConstraints(nTimesteps) +
            ABMConstraints<PredPreyAgent<GRIDSIZE>>::continuityConstraints(nTimesteps);

    std::cout << "Constraints:\n" << fermi << std::endl;

    glp::Problem lp = fermi.toLPProblem();
    std::cout << "LP problem:\n" << lp << std::endl;

    SimplexMCMC  simplex(lp, nullPMF);

    std::cout << "Simplex:\n" << simplex << std::endl;

    int nFractional = 0;
    int nSamples = 1000;
    for(int sample=0; sample < nSamples; ++sample) {
        simplex.randomWalk();
        if(!simplex.solutionIsInteger()) {
            ++nFractional;
            for(int eventId = 1; eventId < simplex.X().size(); ++eventId) {
                if(fabs(simplex.X()[eventId]) > tol)
                    std::cout << simplex.X()[eventId] << " " << Event<PredPreyAgent<GRIDSIZE>>(eventId) << std::endl;
            }
        }
        std::cout << "Sample is: " << nFractional << " " << simplex.X() << std::endl;
    }

//    std::cout << "Simplex:\n" << simplex << std::endl;
    std::cout << "Fractional proportion = " << nFractional * 1.0/nSamples << std::endl;
}




/////////////////////////////////////////////////////////////////////////////////////////
// Creates a low-dimensional measure of a Predator-Prey trajectory in the following way:
//  - calculate the end state of the trajectory
//  - take the total predator and prey occupation of the lower-left corner (rounded down in size)
//    of the grid as two measures and recurse on the upper-right partition until the
//    partition size is smaller than 2x2.
/////////////////////////////////////////////////////////////////////////////////////////
template<int GRIDSIZE>
std::valarray<double> Experiments::Synopsis(const Trajectory<PredPreyAgent<GRIDSIZE>> &trajectory) {
    std::valarray<double> synopsis(floor(log2(GRIDSIZE)) -1);
    ModelState<PredPreyAgent<GRIDSIZE>> endState = trajectory.endState();
    int origin = 0;
    int varid = 0;
    for(int partitionSize = GRIDSIZE / 2; partitionSize > 1; partitionSize /=2) {
        double predOccupation = 0.0;
        double preyOccupation = 0.0;
        for(int x=0; x < partitionSize; ++x) {
            for(int y=0; y < partitionSize; ++y) {
                predOccupation += endState[PredPreyAgent<GRIDSIZE>(origin + x,origin + y,PredPreyAgent<GRIDSIZE>::PREDATOR)];
                preyOccupation += endState[PredPreyAgent<GRIDSIZE>(origin + x,origin + y,PredPreyAgent<GRIDSIZE>::PREY)];
            }
        }
        assert(varid < synopsis.size());
        synopsis[varid++] = predOccupation + preyOccupation;
//        synopsis[varid++] = preyOccupation;
        origin += partitionSize;
    }
    return synopsis;
}

//std::vector<double> Experiments::informationIncrease(int argc, char *argv[]) {
//    if(argc != 9) {
//        std::cout << "Wrong number of arguments. Should be <GRIDSIZE> <nTimestepsPerWindow> <nWindows> <pPredator> <pPrey> <pMakeObservation> <pObserveIfPresent> <nSamplesPerWindow> <nBurnInSamples>" << std::endl;
//        return std::vector<double>();
//    }
//    glp_term_out(GLP_OFF); // turn off GLPK terminal output
//    int gridsize = atoi(argv[1]);
//    int windowSize = atoi(argv[2]);
//    int nWindows = atoi(argv[3]);
//    double pPredator = atof(argv[4]);         // Poisson prob of predator in each gridsquare at t=0
//    double pPrey = atof(argv[5]);             // Poisson prob of prey in each gridsquare at t=0
//    double pMakeObservation = atof(argv[6]);  // prob of making an observation of each gridsquare at each timestep
//    double pObserveIfPresent = atof(argv[7]); // prob of detecting an agent given that it is present
//    int nSamplesPerWindow = atoi(argv[8]);
//    int nBurnInSamples = atoi(argv[9]);
//    return informationIncrease(gridsize, windowSize, nWindows, pPredator, pPrey,
//                               pMakeObservation, pObserveIfPresent, nSamplesPerWindow, nBurnInSamples);
//}


template<int GRIDSIZE>
std::vector<double> Experiments::informationIncrease(
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

//    std::vector<boost::math::binomial_distribution<double>> startWeights(PredPreyAgent::domainSize());
//    for(int agentId=0; agentId < PredPreyAgent::domainSize(); ++agentId) {
//        startWeights[agentId] = boost::math::binomial_distribution<double>(
//                1,(PredPreyAgent(agentId).type() == PredPreyAgent::PREDATOR?pPredator:pPrey)
//                );
//    }

    BernoulliModelState<PredPreyAgent<GRIDSIZE>> startStateDist([pPredator,pPrey](PredPreyAgent<GRIDSIZE> agent) {
        return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
    });

    AssimilationWindow<PredPreyAgent<GRIDSIZE>> window(windowSize, startStateDist, pMakeObservation, pObserveIfPresent);

    MCMCSampler sampler(window.posterior);
//    solver.solve(nBurnInSamples, nSamplesPerWindow);


    return std::vector<double>(); //assimilation.calculateInformationGain();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


//Gnuplot &Experiments::plotAgents(Gnuplot &gp, const ModelState<PredPreyAgent> &state) {
//    std::vector<std::tuple<double,double,double>> pointData;
//    for(int x=0; x<PredPreyAgent::GRIDSIZE; ++x) {
//        for(int y=0; y<PredPreyAgent::GRIDSIZE; ++y) {
//            int colour = 2*(state[PredPreyAgent(x,y,PredPreyAgent::PREDATOR)]>0.0)
//                         + (state[PredPreyAgent(x,y,PredPreyAgent::PREY)]>0.0);
//            if(colour != 0)
//                pointData.emplace_back(x, y, colour);
//        }
//    }
//
//    gp << "set linetype 1 lc 'red'\n";
//    gp << "set linetype 2 lc 'blue'\n";
//    gp << "plot [-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE-0.5 << "] ";
//    gp << "'-' with points pointtype 5 pointsize 0.5 lc variable\n";
//    gp.send1d(pointData);
//    return gp;
//}
//
