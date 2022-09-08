//
// Created by daniel on 06/05/2021.
//

#include "Experiments.h"
#include "BernoulliStartState.h"
#include "Likelihood.h"
#include "SparseBasisSampler.h"
#include "agents/BinomialAgent.h"
#include "RejectionSampler.h"
#include "diagnostics/Dataflow.h"
#include "diagnostics/AgentDataflow.h"
#include "agents/CatMouseAgent.h"
#include "agents/PredPreyAgent.h"
#include "Prior.h"
#include "ExactSolver.h"
#include "PoissonStartState.h"

using namespace dataflow;

void Experiments::BinomialAgentSingleObservation() {
    constexpr int nTimesteps = 2;
    constexpr int GRIDSIZE = nTimesteps+1;
    constexpr int nSamples = 100000;
    constexpr int nBurnin = 100;
    constexpr int nRejectionSamples = 50000;
    constexpr double kappa = 1.0;

//    doSingleObservationExperiment(nTimesteps, nBurnin, nSamples, nRejectionSamples, kappa,
//                                  BernoulliStartState<BinomialAgent<GRIDSIZE>>({1.0, 0.1, 0.0}),
//                                  AgentStateObservation(State<BinomialAgent<GRIDSIZE>>(1, 0),1,0.9));
    std::cout << "Exact state = " << 0.5+0.1*0.25 << " " << 0.5 + 0.1*0.25 << " " << 0.1*0.5 << std::endl;
}


void Experiments::CatMouseSingleObservation() {
    constexpr int nTimesteps = 2;
    constexpr int nBurnin = 100;
    constexpr int nSamples = 200000;
    constexpr int nRejectionSamples = 200000;
    constexpr double kappa = 1.75;

//    doSingleObservationExperiment(nTimesteps, nBurnin, nSamples, nRejectionSamples, kappa,
//            BernoulliStartState<CatMouseAgent>({0.9, 0.1, 0.1, 0.9}),
//            AgentStateObservation(State(1,CatMouseAgent(CatMouseAgent::CAT, CatMouseAgent::LEFT)),1,1.0));
}


void Experiments::PredPreySingleObservation() {
    constexpr int GRIDSIZE = 3;
    constexpr int nTimesteps = 2;
    constexpr double pPredator = 0.1;//0.08;          // Bernoulli prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Bernoulli prob of prey in each gridsquare at t=0
    constexpr int nSamples = 500000; //250000;
    constexpr int nBurnin = 10000;
    constexpr int nRejectionSamples = 100000;
    constexpr double kappa = 3.0;

//    doSingleObservationExperiment(nTimesteps, nBurnin, nSamples, nRejectionSamples, kappa,
//            PoissonStartState<PredPreyAgent<GRIDSIZE>>([](PredPreyAgent<GRIDSIZE> agent) {
//                return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
//            }),
//            AgentStateObservation(State(1,PredPreyAgent<GRIDSIZE>(1, 1, PredPreyAgent<GRIDSIZE>::PREY)),1,1.0)
//    );

//    Plotter gp;
//    gp.plot(window.realTrajectory.endState(), sampleStats.means());
//
//    Plotter gp2;
//    gp2.plot(window.realTrajectory.endState(), rejectionSampleStats.means());

}


void Experiments::CatMouseAssimilation() {
    constexpr int nTimesteps = 2;
    constexpr double pMakeObservation = 0.2;
    constexpr double pObserveIfPresent = 1.0;
    constexpr int nSamples = 200000;
    constexpr int nBurnin = 100;
    constexpr double kappa = 1.25;

    Prior<CatMouseAgent> prior(nTimesteps, PoissonStartState<CatMouseAgent>({0.5, 0.5, 0.3, 0.3}));
    std::cout << "Prior support is\n" << prior << std::endl;

    Trajectory<CatMouseAgent> realTrajectory = prior.nextSample();

    Likelihood<CatMouseAgent> likelihood(realTrajectory, pMakeObservation, pObserveIfPresent);
    std::cout << "Likelihood support is\n" << likelihood << std::endl;

    auto posterior = likelihood * prior;
//
//    SparseBasisSampler sampler(posterior,kappa);
//    std::cout << "Constructed basis\n" << sampler << std::endl;
//    for(int i=0; i < nBurnin; ++i) {
//        const std::vector<ABM::occupation_type> &sample = sampler();
//        assert(posterior.constraints.isValidSolution(sample));
//    }
//    std::map<std::vector<ABM::occupation_type>,double> pmf;
//    for(int s=0; s<nSamples; ++s) pmf[Trajectory<CatMouseAgent>(sampler(),nTimesteps)] += 1.0/nSamples;
//    std::cout << sampler.stats;
//
//    RejectionSampler<Trajectory<CatMouseAgent>> rejectionSampler(prior,likelihood);
//    std::map<std::vector<ABM::occupation_type>,double> rejectionPMF;
//    for(int s=0; s<nSamples; ++s) {
//        Trajectory<CatMouseAgent> trajectory(rejectionSampler(),nTimesteps);
////        std::cout << "Got trajectory " << trajectory << std::endl;
//        rejectionPMF[trajectory] += 1.0/nSamples;
//    }
//
//    ExactSolver<CatMouseAgent> exactSolver(posterior);
//    std::cout << "Comparison of exact / rejection / mcmc solutions:" << std::endl;
//    std::fixed(std::cout).precision(5);
//    for(auto &[trajectory,prob] : exactSolver.pmf) {
//        std::cout << trajectory << "  p(ex)= " << prob << "  p(rej)= " << rejectionPMF[trajectory] << "  p(mcmc)= " << pmf[trajectory] << "  err= " << prob - pmf[trajectory] << std::endl;
//    }
}




template<class AGENT>
void Experiments::doSingleObservationExperiment(int nTimesteps, int nBurnin, int nSamples, int nRejectionSamples,
                                                double kappa,
                                                const StartStateDistribution<AGENT> &startState,
                                                const NoisyAgentStateObservation<AGENT> &observation) {
//    Prior<AGENT> prior(nTimesteps, startState);
//    std::cout << "Prior support is\n" << prior << std::endl;
//
//    Likelihood<AGENT> likelihood(observation);
//    std::cout << "Likelihood support is\n" << likelihood << std::endl;
//
//    auto posterior = likelihood * prior;
//    std::cout << "Posterior support is\n" << posterior << std::endl;
//
//    SparseBasisSampler sampler(posterior, kappa);
//    std::cout << "Constructed basis\n" << sampler << std::endl;
//
//    for(int i=0; i < nBurnin; ++i) {
//        const std::vector<ABM::occupation_type> &sample = sampler();
////        std::cout << "Got sample " << sample << "  " << ModelState<AGENT>(sample, nTimesteps, nTimesteps) << std::endl;
//        assert(posterior.constraints.isValidSolution(sample));
//    }
//
//    std::fixed(std::cout).precision(5);
//    ModelState<AGENT> aggregateState;
//    sampler >>= Take(nSamples) >>= TrajectoryToModelState<AGENT>(nTimesteps, nTimesteps) >>= Sum(aggregateState);
//    std::cout << sampler.stats << std::endl;
//    std::cout << "     MCMC state = " << aggregateState / nSamples << std::endl;
//
//    RejectionSampler<Trajectory<AGENT>> rejectionSampler(prior,likelihood);
//
//    ModelState<AGENT> rejectionAggregateState;
//    rejectionSampler >>= Take(nRejectionSamples) >>= TrajectoryToModelState<AGENT>(nTimesteps, nTimesteps) >>= Sum(rejectionAggregateState);
//    std::cout << "Rejection state = " << rejectionAggregateState / nRejectionSamples << std::endl;
//    std::valarray<double> err = (aggregateState / nSamples) - (rejectionAggregateState / nRejectionSamples);
//    double rms = sqrt((err*err).sum()/err.size());
//    std::cout << "      RMS Error = " << rms << std::endl;
}
