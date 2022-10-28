//
// Created by daniel on 06/05/2021.
//

#include "Experiments.h"
#include "BernoulliStartState.h"
#include "ABMLikelihood.h"
#include "agents/BinomialAgent.h"
#include "RejectionSampler.h"
#include "diagnostics/Dataflow.h"
#include "agents/CatMouseAgent.h"
#include "agents/PredPreyAgent.h"
#include "ABMPrior.h"
#include "ExactSolver.h"
#include "PoissonStartState.h"
#include "FactorisedDistributionSampler.h"
#include "FixedPopulationStartState.h"
#include "ExtendedTrajectory.h"
#include "PredPreyTrajectory.h"
#include "Plotter.h"
#include "ABMPriorSampler.h"
#include "diagnostics/MultiChainStats.h"
#include "include/asyncvector.h"
#include "agents/ReducedPredPreyAgent.h"
#include "TableauEntropyMaximiser.h"

using namespace dataflow;

void Experiments::BinomialAgentSingleObservation() {
    constexpr int nTimesteps = 2;
    constexpr int GRIDSIZE = nTimesteps+1;
    constexpr int nSamples = 100000;
    constexpr int nBurnin = 100;
    constexpr int nRejectionSamples = 100000;
    constexpr double kappa = 1.75;

    double p0 = 1.0;
    double p1 = 0.1;
    double pObserveIfPresent = 0.9;

    typedef BinomialAgent<GRIDSIZE> agent_type;
    typedef ExtendedTrajectory<agent_type,nTimesteps> trajectory_type;

    BernoulliStartState<agent_type> startState({p0, p1, 0.0}, kappa);
    ABMLikelihood<trajectory_type> likelihood(State<agent_type>(1,0), 1, pObserveIfPresent, kappa);
    ABMPosterior posterior(startState, likelihood);

    std::cout << posterior.constraints << std::endl;

    std::vector<double> entropies(trajectory_type::dimension);
    for(double &EFi: entropies) EFi = Random::nextDouble();

    TableauEntropyMaximiser<int> tableau(entropies, posterior.constraints, 0.0001);
    std::cout << tableau << std::endl;

    tableau.factorise();

    std::cout << tableau << std::endl;

//    doValidationExperiment(posterior, nBurnin, nSamples, nRejectionSamples);
//    std::cout << "Exact state = " << 0.5*p0+ 0.25*p1 << " " << 0.5*p0 + 0.25*p1 << " " << 0.5*p1 << std::endl;
}


void Experiments::CatMouseSingleObservation() {
    constexpr int nTimesteps = 2;
    constexpr int nBurnin = 500;
    constexpr int nSamples = 500000;
    constexpr int nRejectionSamples = 100000;
    constexpr double kappa = 1.75;
    constexpr double pObserveIfPresent = 1.0;

//    doValidationExperiment(nTimesteps, nBurnin, nSamples, nRejectionSamples, kappa,
//            BernoulliStartState<Trajectory<CatMouseAgent>>({0.9, 0.1, 0.1, 0.9}),
//            Likelihood<Trajectory<CatMouseAgent>>(State(1,CatMouseAgent(CatMouseAgent::CAT, CatMouseAgent::LEFT)),1,1.0));

    typedef CatMouseAgent agent_type;
    typedef ExtendedTrajectory<agent_type,nTimesteps> trajectory_type;

    BernoulliStartState<agent_type> startState({0.9, 0.1, 0.1, 0.9}, kappa);
    ABMLikelihood<trajectory_type> likelihood(State(1, CatMouseAgent(CatMouseAgent::CAT, CatMouseAgent::LEFT)), 1, pObserveIfPresent, kappa);
    ABMPosterior posterior(startState, likelihood);

    doValidationExperiment(posterior, nBurnin, nSamples, nRejectionSamples);
}


void Experiments::PredPreySingleObservation() {
    constexpr int GRIDSIZE = 3;
    constexpr int nTimesteps = 2;
    constexpr double pPredator = 0.1;//0.08;          // Bernoulli prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Bernoulli prob of prey in each gridsquare at t=0
    constexpr int nSamples = 100000; //250000;
    constexpr int nBurnin = nSamples/5;
    constexpr int nRejectionSamples = nSamples/10;
    constexpr double kappa = 3.5;// 5.75;

    typedef PredPreyAgent<GRIDSIZE> agent_type;
    typedef PredPreyTrajectory<GRIDSIZE,nTimesteps> trajectory_type;

    PoissonStartState<agent_type> startState([](agent_type agent) {
        return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
    }, kappa);
    ABMLikelihood<trajectory_type> likelihood(State(1, PredPreyAgent<GRIDSIZE>(1, 1, PredPreyAgent<GRIDSIZE>::PREY)), 1, 1.0, kappa);
    ABMPosterior posterior(startState, likelihood);

    doValidationExperiment(posterior, nBurnin, nSamples, nRejectionSamples);

}


void Experiments::PredPreyPriorTest() {
    constexpr int GRIDSIZE = 8;
    constexpr int TIMESTEPS = 4;
    constexpr double pPredator = 0.05;//0.08;          // Bernoulli prob of predator in each gridsquare at t=0
    constexpr double pPrey = pPredator;    // Bernoulli prob of prey in each gridsquare at t=0
    constexpr int nSamples = 100000; //250000;
    constexpr double kappa = 6.0;//4.75;//3.5;

    typedef PredPreyAgent<GRIDSIZE> agent_type;
    typedef PredPreyTrajectory<GRIDSIZE,TIMESTEPS> trajectory_type;

    PoissonStartState<agent_type> startState([](agent_type agent) {
        return agent.type() == agent_type::PREDATOR?pPredator:pPrey;
    }, kappa);

    ABMPrior prior = makeABMPrior<trajectory_type>(startState);

    std::cout << prior;

//    FactorisedDistributionSampler sampler(prior);
//        std::cout << "Burning-in..." << std::endl;
//        for(int s = 1; s<1000; ++s) {
//            sampler(); // burn-in
//        }
//        std::cout << "Sampling..." << std::endl;
//        MultiChainStats multiChainStats(1000000,sampler);
//    MultiChainStats result(nSamples/2, prior);
//    std::cout << result;


//    MultiChainStats multiChainStats(generateVectorAsync(1, [prior]() {
//        return MultiChainStats(nSamples/2, prior);
//    }));

    MultiChainStats multiChainStats(generateVectorAsync(1, [&prior, nSamples]() {
        FactorisedDistributionSampler sampler(prior);
        std::cout << sampler.basisVectors << std::endl;
        std::cout << "Burning-in..." << std::endl;
        for(int s = 1; s<nSamples/5; ++s) sampler(); // burn-in
        std::cout << sampler.stats << std::endl;
        std::cout << "Sampling..." << std::endl;
        return MultiChainStats(nSamples, sampler);
    }));
//
    Plotter plotter;
    plotter.plot<agent_type>(multiChainStats);

}


void Experiments::ReducedPredPreyPriorTest() {
    constexpr int GRIDSIZE = 8;
    constexpr int TIMESTEPS = 4;
    constexpr double pPredator = 0.05;//0.08;          // Bernoulli prob of predator in each gridsquare at t=0
    constexpr double pPrey = pPredator;    // Bernoulli prob of prey in each gridsquare at t=0
    constexpr int nSamples = 500000; //250000;
    constexpr double kappa = 6.0;//4.75;//3.5;

    typedef ReducedPredPreyAgent<GRIDSIZE> agent_type;
    typedef ExtendedTrajectory<agent_type,TIMESTEPS> trajectory_type;

    PoissonStartState<agent_type> startState([](agent_type agent) {
        return agent.type() == agent_type::PREDATOR?pPredator:pPrey;
    }, kappa);

    ABMPrior prior = makeABMPrior<trajectory_type>(startState);

    std::cout << prior;

//    FactorisedDistributionSampler sampler(prior);
//        std::cout << "Burning-in..." << std::endl;
//        for(int s = 1; s<1000; ++s) {
//            sampler(); // burn-in
//        }
//        std::cout << "Sampling..." << std::endl;
//        MultiChainStats multiChainStats(1000000,sampler);
//    MultiChainStats result(nSamples/2, prior);
//    std::cout << result;


    MultiChainStats multiChainStats(generateVectorAsync(4, [&prior, nSamples]() {
        FactorisedDistributionSampler sampler(prior);
        std::cout << "Burning-in..." << std::endl;
        for(int s = 1; s<nSamples/5; ++s) sampler(); // burn-in
        std::cout << sampler.stats << std::endl;
        std::cout << "Sampling..." << std::endl;
        return MultiChainStats(nSamples, sampler);
    }));
//
    Plotter plotter;
    plotter.plot<agent_type>(multiChainStats);

}


//void Experiments::CatMouseAssimilation() {
//    constexpr int nTimesteps = 2;
//    constexpr double pMakeObservation = 0.2;
//    constexpr double pObserveIfPresent = 1.0;
//    constexpr int nSamples = 1000000;
//    constexpr int nBurnin = 100;
//    constexpr double kappa = 3.0;
//
//
//
//    Prior<Trajectory<CatMouseAgent>> prior(nTimesteps, FixedPopulationStartState<Trajectory<CatMouseAgent>>({0, 0, 1, 1},{1,1}));
//    std::cout << "Prior is\n" << prior << std::endl;
//
//    Trajectory<CatMouseAgent> realTrajectory = prior.sampler()();
//    std::cout << "Real trajectory " << realTrajectory << std::endl;
//
//    Likelihood<Trajectory<CatMouseAgent>> likelihood(realTrajectory, pMakeObservation, pObserveIfPresent);
//    std::cout << "Likelihood is\n" << likelihood << std::endl;
//
//    auto posterior = likelihood * prior;
//    std::cout << "Posterior is\n" << posterior << std::endl;
//
//    Trajectory<CatMouseAgent> initialSample = prior.sampler()();
//    std::cout << "Initial sample " << initialSample << std::endl;
//
//    debug(posterior.sanityCheck(initialSample));
//    assert(posterior.isFeasible(realTrajectory));
//
//    ConstrainedFactorisedSampler<Trajectory<CatMouseAgent>> sampler(posterior, initialSample);
//    debug(sampler.sanityCheck());
//
//    std::cout << "Burning in..." << std::endl;
//    for(int i=0; i < nBurnin; ++i) {
//        const Trajectory<CatMouseAgent> &sample = sampler.nextSample();
//        assert(posterior.constraints.isValidSolution(sample));
////        std::cout << "Infeasibility = " << modelStateSampler.currentInfeasibility << std::endl;
//    }
//    std::cout << "Sampling..." << std::endl;
//    std::map<std::vector<ABM::occupation_type>,double> pmf;
//    int reportInterval = nSamples/10;
//    for(int s=0; s<nSamples; ++s) {
//        if(s%reportInterval == 0) std::cout << s*100.0/nSamples << "% complete" << std::endl;
//        pmf[sampler.nextSample()] += 1.0/nSamples;
//    }
//    std::cout << sampler.stats;
//
//    std::cout << "Validating against rejection sampler..." << std::endl;
//    RejectionSampler<Trajectory<CatMouseAgent>> rejectionSampler(prior,likelihood);
//    std::map<std::vector<ABM::occupation_type>,double> rejectionPMF;
//    for(int s=0; s<nSamples; ++s) {
//        Trajectory<CatMouseAgent> trajectory(rejectionSampler(),nTimesteps);
////        std::cout << "Got trajectory " << trajectory << std::endl;
//        rejectionPMF[trajectory] += 1.0/nSamples;
//    }
//
//    std::cout << "Validating against exact solution..." << std::endl;
//    ExactSolver<Trajectory<CatMouseAgent>> exactSolver(posterior, Trajectory<CatMouseAgent>(nTimesteps));
//    std::cout << "Comparison of exact / rejection / mcmc solutions:" << std::endl;
//    std::fixed(std::cout).precision(5);
//    for(auto &[trajectory,prob] : exactSolver.pmf) {
//        std::cout << trajectory << "  p_exact = " << prob << "  p_reject = " << rejectionPMF[trajectory] << "  p_mcmc = " << pmf[trajectory] << "  p_exact - p_mcmc = " << prob - pmf[trajectory] << std::endl;
//    }
//}




template<class DOMAIN, class STARTSTATE>
void Experiments::doValidationExperiment(const ABMPosterior<DOMAIN,STARTSTATE> &posterior, int nBurnin, int nSamples, int nRejectionSamples) {
//    std::cout << "Prior is\n" << posterior.prior << std::endl;
//
//    std::cout << "Likelihood is\n" << posterior.likelihood << std::endl;
//
    std::cout << "Posterior is\n" << posterior << std::endl;

    FactorisedDistributionSampler sampler(posterior);

    std::cout << "Burning-in..." << std::endl;
    for(int i=0; i < nBurnin; ++i) {
        const DOMAIN &sample = sampler.nextSample();
//        std::cout << "Got sample " << sample << "  " << sample.endState() << std::endl;
        assert(posterior.constraints.isValidSolution(sample));
    }

    std::cout << "Sampling..." << std::endl;
    std::fixed(std::cout).precision(5);
    std::valarray<double> aggregateState(0.0,DOMAIN::agent_type::domainSize);
    auto startTime = std::chrono::steady_clock::now();
    sampler >>= &DOMAIN::endState >>= Sum(nSamples,aggregateState);
    auto endTime = std::chrono::steady_clock::now();
    aggregateState /= nSamples;
    std::cout << "Sample time = " << endTime - startTime << std::endl;
    std::cout << sampler.stats << std::endl;
    std::cout << "     MCMC state = " << aggregateState << std::endl;

    RejectionSampler<DOMAIN> rejectionSampler(posterior.prior, posterior.likelihood);
    std::valarray<double> rejectionAggregateState(0.0,DOMAIN::agent_type::domainSize);
    rejectionSampler >>= &DOMAIN::endState >>= Sum(nRejectionSamples, rejectionAggregateState);
    rejectionAggregateState /= nRejectionSamples;
    std::cout << "Rejection state = " << rejectionAggregateState << std::endl;

    std::valarray<double> err = aggregateState - rejectionAggregateState;
    double rms = sqrt((err*err).sum()/err.size());
    std::cout << " RMS difference = " << rms << std::endl;
}

//void Experiments::Animation() {
//    const int GRIDSIZE = 32;
//    const int nTimesteps = 16;
//
////    std::cout << exp(PredPreyAgentBase::lpPredBirthGivenPrey) << std::endl;
////    std::cout << exp(PredPreyAgentBase::lpPredDeathGivenPrey) << std::endl;
////    std::cout << exp(PredPreyAgentBase::lpPredBirthGivenNoPrey) << std::endl;
////    std::cout << exp(PredPreyAgentBase::lpPredDeathGivenNoPrey) << std::endl;
////    std::cout << std::endl;
////    std::cout << exp(PredPreyAgentBase::lpPreyBirthGivenPred) << std::endl;
////    std::cout << exp(PredPreyAgentBase::lpPreyDeathGivenPred) << std::endl;
////    std::cout << exp(PredPreyAgentBase::lpPreyBirthGivenNoPred) << std::endl;
////    std::cout << exp(PredPreyAgentBase::lpPreyDeathGivenNoPred) << std::endl;
////    std::cout << std::endl;
////    std::cout << exp(PredPreyAgentBase::lpMove) << std::endl;
//
//    PoissonStartState<Trajectory<PredPreyAgent<GRIDSIZE>>> startState([](PredPreyAgent<GRIDSIZE> agent) {
//        return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?PredPreyAgentBase::pPred:PredPreyAgentBase::pPrey;
//    });
//
//    Prior<Trajectory<PredPreyAgent<GRIDSIZE>>> prior(nTimesteps, startState);
//
//    std::cout << "Sampling" << std::endl;
//    Trajectory trajectory = prior.sampler()();
//    std::cout << "Got sample" << std::endl;
//
//    Plotter plotter;
//
//    std::cout << "Plotting..." << std::endl;
//    plotter.animate(trajectory, 15);
//    std::cout << "Done" << std::endl;
////    std::cout << "Prior is\n" << prior << std::endl;
////
////    Trajectory<CatMouseAgent> realTrajectory = prior.modelStateSampler()();
////    std::cout << "Real trajectory " << realTrajectory << std::endl;
//
//}