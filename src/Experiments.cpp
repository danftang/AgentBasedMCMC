//
// Created by daniel on 06/05/2021.
//

#include "Experiments.h"
#include "BernoulliStartState.h"
#include "Likelihood.h"
#include "agents/BinomialAgent.h"
#include "RejectionSampler.h"
#include "diagnostics/Dataflow.h"
#include "agents/CatMouseAgent.h"
#include "agents/PredPreyAgent.h"
#include "ABMPrior.h"
#include "ExactSolver.h"
#include "PoissonStartState.h"
#include "ConstrainedFactorisedSampler.h"
#include "FixedPopulationStartState.h"
#include "ExtendedTrajectory.h"
#include "PredPreyTrajectory.h"
#include "Plotter.h"
#include "ABMPriorSampler.h"

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

    doValidationExperiment(nBurnin, nSamples, nRejectionSamples, kappa,
                           BernoulliStartState<ExtendedTrajectory<BinomialAgent<GRIDSIZE>,nTimesteps>>({p0, p1, 0.0}),
                           Likelihood<ExtendedTrajectory<BinomialAgent<GRIDSIZE>,nTimesteps>>(State<BinomialAgent<GRIDSIZE>>(1, 0), 1, pObserveIfPresent));
//    doValidationExperiment(nTimesteps, nBurnin, nSamples, nRejectionSamples, kappa,
//                           BernoulliStartState<ExtendedTrajectory<BinomialAgent<GRIDSIZE>>>({p0, p1, 0.0}),
//                           Likelihood<ExtendedTrajectory<BinomialAgent<GRIDSIZE>>>(State<BinomialAgent<GRIDSIZE>>(1, 0), 1, pObserveIfPresent));
    std::cout << "Exact state = " << 0.5*p0+ 0.25*p1 << " " << 0.5*p0 + 0.25*p1 << " " << 0.5*p1 << std::endl;
}


void Experiments::CatMouseSingleObservation() {
    constexpr int nTimesteps = 2;
    constexpr int nBurnin = 100;
    constexpr int nSamples = 100000;
    constexpr int nRejectionSamples = 100000;
    constexpr double kappa = 2.25;//1.75;

//    doValidationExperiment(nTimesteps, nBurnin, nSamples, nRejectionSamples, kappa,
//            BernoulliStartState<Trajectory<CatMouseAgent>>({0.9, 0.1, 0.1, 0.9}),
//            Likelihood<Trajectory<CatMouseAgent>>(State(1,CatMouseAgent(CatMouseAgent::CAT, CatMouseAgent::LEFT)),1,1.0));
    doValidationExperiment(nBurnin, nSamples, nRejectionSamples, kappa,
                           BernoulliStartState<ExtendedTrajectory<CatMouseAgent,nTimesteps>>({0.9, 0.1, 0.1, 0.9}),
                           Likelihood<ExtendedTrajectory<CatMouseAgent,nTimesteps>>(State(1, CatMouseAgent(CatMouseAgent::CAT, CatMouseAgent::LEFT)), 1, 1.0));
}


void Experiments::PredPreySingleObservation() {
    constexpr int GRIDSIZE = 3;
    constexpr int nTimesteps = 2;
    constexpr double pPredator = 0.1;//0.08;          // Bernoulli prob of predator in each gridsquare at t=0
    constexpr double pPrey = 2.0*pPredator;    // Bernoulli prob of prey in each gridsquare at t=0
    constexpr int nSamples = 200000; //250000;
    constexpr int nBurnin = 1000;
    constexpr int nRejectionSamples = 100000;
    constexpr double kappa = 5.0;// 5.75;

//    doValidationExperiment(nBurnin, nSamples, nRejectionSamples, kappa,
//                           PoissonStartState<ExtendedTrajectory2<PredPreyAgent<GRIDSIZE>,nTimesteps>>([](PredPreyAgent<GRIDSIZE> agent) {
//                               return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
//                           }),
//                           Likelihood<ExtendedTrajectory2<PredPreyAgent<GRIDSIZE>,nTimesteps>>(State(1,PredPreyAgent<GRIDSIZE>(1, 1, PredPreyAgent<GRIDSIZE>::PREY)),1,1.0)
//    );

    doValidationExperiment(nBurnin, nSamples, nRejectionSamples, kappa,
                           PoissonStartState<PredPreyTrajectory<GRIDSIZE,nTimesteps>>([](PredPreyAgent<GRIDSIZE> agent) {
                               return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
                           }),
                           Likelihood<PredPreyTrajectory<GRIDSIZE,nTimesteps>>(State(1,PredPreyAgent<GRIDSIZE>(1, 1, PredPreyAgent<GRIDSIZE>::PREY)),1,1.0)
    );

//    Plotter gp;
//    gp.plot(window.realTrajectory.endState(), sampleStats.means());
//
//    Plotter gp2;
//    gp2.plot(window.realTrajectory.endState(), rejectionSampleStats.means());

}


//void Experiments::CatMouseAssimilation() {
//    constexpr int nTimesteps = 2;
//    constexpr double pMakeObservation = 0.2;
//    constexpr double pObserveIfPresent = 1.0;
//    constexpr int nSamples = 1000000;
//    constexpr int nBurnin = 100;
//    constexpr double kappa = 3.0;
//
//    ABM::kappa = kappa;
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
//
//


template<class DOMAIN, class STARTSTATE>
void Experiments::doValidationExperiment(int nBurnin, int nSamples, int nRejectionSamples,
                                         double kappa,
                                         const STARTSTATE &startStateDistribution,
                                         const ConstrainedFactorisedDistribution<DOMAIN> &likelihood) {
    ABM::kappa = kappa;
    ABMPrior<DOMAIN> prior(startStateDistribution);
    std::cout << "Prior is\n" << prior << std::endl;

    std::cout << "Likelihood is\n" << likelihood << std::endl;

    auto posterior = likelihood * prior;
    std::cout << "Posterior is\n" << posterior << std::endl;


    ConstrainedFactorisedSampler sampler(posterior);

    std::cout << "Burning-in..." << std::endl;
    for(int i=0; i < nBurnin; ++i) {
        const DOMAIN &sample = sampler.nextSample();
//        std::cout << "Got sample " << sample << "  " << sample.endState() << std::endl;
        assert(posterior.constraints.isValidSolution(sample));
    }

    std::cout << "Sampling..." << std::endl;
    std::fixed(std::cout).precision(5);
    ModelState<typename DOMAIN::agent_type> aggregateState;
    auto startTime = std::chrono::steady_clock::now();
    sampler >>= &DOMAIN::endState >>= Sum(nSamples,aggregateState);
    auto endTime = std::chrono::steady_clock::now();
    std::cout << "Sample time = " << endTime - startTime << std::endl;
    std::cout << sampler.stats << std::endl;
    std::cout << "     MCMC state = " << aggregateState / nSamples << std::endl;

    RejectionSampler<DOMAIN> rejectionSampler(startStateDistribution.priorSampler(), likelihood);
    ModelState<typename DOMAIN::agent_type> rejectionAggregateState;
    rejectionSampler >>= &DOMAIN::endState >>= Sum(nRejectionSamples, rejectionAggregateState);
    std::cout << "Rejection state = " << rejectionAggregateState / nRejectionSamples << std::endl;

    std::valarray<double> err = (aggregateState / nSamples) - (rejectionAggregateState / nRejectionSamples);
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