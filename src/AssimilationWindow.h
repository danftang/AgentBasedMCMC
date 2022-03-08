//
// Created by daniel on 20/07/2021.
//

#ifndef GLPKTEST_ASSIMILATIONWINDOW_H
#define GLPKTEST_ASSIMILATIONWINDOW_H

#include <vector>
#include "debug.h"
#include "ConvexPMF.h"
#include "TrajectoryLikelihoodPMF.h"
#include "TrajectoryPMF.h"
#include "TrajectorySampler.h"
#include "BinomialDistribution.h"
#include "gnuplot-iostream/gnuplot-iostream.h"
//#include "agents/PredPreyAgent.h"
#include "RejectionSampler.h"

// For testing puroposes only
template<typename AGENT>
class AssimilationWindow {
public:
    ConvexPMF<Trajectory<AGENT>>        priorPMF;
    std::function<Trajectory<AGENT>()>  priorSampler;
    Trajectory<AGENT>                   realTrajectory;
    ConvexPMF<Trajectory<AGENT>>        likelihoodPMF;
    ConvexPMF<Trajectory<AGENT>>        posterior;


    AssimilationWindow(
            const Distribution<ModelState<AGENT>> &priorStartState,
            Trajectory<AGENT> realTraj,
            double pMakeObservation,
            double pObserveIfPresent)
            :
            priorPMF(realTraj.nTimesteps(), priorStartState.PMF()),
            priorSampler(Trajectory<AGENT>::priorSampler(realTraj.nTimesteps(), priorStartState.sampler())),
            realTrajectory(std::move(realTraj)),
            likelihoodPMF(realTrajectory, pMakeObservation, pObserveIfPresent),
            posterior(likelihoodPMF * priorPMF)
    { }


    AssimilationWindow(
            int nTimesteps,
            ConvexPMF<ModelState<AGENT>> priorStartStatePMF,
            std::function<ModelState<AGENT>()> priorStartStateSampler,
            double pMakeObservation,
            double pObserveIfPresent)
    :
    priorPMF(nTimesteps, std::move(priorStartStatePMF)),
    priorSampler(Trajectory<AGENT>::priorSampler(nTimesteps, std::move(priorStartStateSampler))),
    realTrajectory(priorSampler()),
    likelihoodPMF(realTrajectory, pMakeObservation, pObserveIfPresent),
    posterior(likelihoodPMF * priorPMF)
    {
//        debug(std::cout << "Created window with real trajectory " << realTrajectory << std::endl);
    }

    AssimilationWindow(int nTimesteps, const Distribution<ModelState<AGENT>> &priorStartState, double pMakeObservation, double pObserveIfPresent)
    : AssimilationWindow(nTimesteps, priorStartState.PMF(), priorStartState.sampler(), pMakeObservation, pObserveIfPresent)
    {
    }


    AssimilationWindow(int nTimesteps, const Distribution<ModelState<AGENT>> &priorStartState, const AgentStateObservation<AGENT> &observation)
    :
    priorPMF(nTimesteps, priorStartState.PMF()),
    priorSampler(Trajectory<AGENT>::priorSampler(nTimesteps, priorStartState.sampler())),
    likelihoodPMF(nTimesteps, observation),
    realTrajectory(0),
    posterior(likelihoodPMF * priorPMF)
    {
        RejectionSampler<Trajectory<AGENT>> rejectionSampler(priorSampler, likelihoodPMF);
        realTrajectory = rejectionSampler();
    }

    AssimilationWindow(int nTimesteps, const Distribution<ModelState<AGENT>> &priorStartState)
    :
    priorPMF(nTimesteps, priorStartState.PMF()),
    priorSampler(Trajectory<AGENT>::priorSampler(nTimesteps, priorStartState.sampler())),
    realTrajectory(Trajectory<AGENT>(nTimesteps,priorSampler())),
    likelihoodPMF(ConvexPMF<Trajectory<AGENT>>::uniformTrajectoryDistribution(nTimesteps)),
    posterior(priorPMF)
    {
    }

    int dimension() { return realTrajectory.dimension(); }
//    void doAnalysis(int nSamples, int nBurnInSamples) {
//        SimplexMCMC sampler = SimplexMCMC(posterior, priorSampler());
//
//        std::cout << "Starting burn-in" << std::endl;
//        for(int burnIn=0; burnIn<nBurnInSamples; ++burnIn) {
//            sampler.nextSample();
//        }
//
////        std::cout << "Sampler:\n" << sampler << std::endl;
//        std::cout << "Initial exactEndState: " << sampler.X() << std::endl;
//        for(int s=0; s<nSamples; ++s) {
//            Trajectory<AGENT> sample(sampler.nextSample());
//            //            std::cout << "Sampler:\n" << sampler << std::endl;
//            //            std::cout << "Sample: " << sample << std::endl;
//            assert(posterior.convexSupport.isValidSolution(sample));
//            analysis += sample.endState();
//        }
//        std::cout << "Feasible stats:\n" << sampler.feasibleStatistics << std::endl;
//        std::cout << "Infeasible stats:\n" << sampler.infeasibleStatistics << std::endl;
//        std::cout << "Infeasible proportion = " << sampler.infeasibleStatistics.nSamples*1.0/(sampler.feasibleStatistics.nSamples + sampler.infeasibleStatistics.nSamples) << std::endl;
//
////        std::cout << analysis << std::endl;
//    }


//////////// Use TrajectorySampler instead!
//    BinomialDistribution priorEndState(int nSamples) const {
//        MCMCStatistics endStats(AGENT::domainSize());
//        for(int s=0; s<nSamples; ++s) {
//            Trajectory<AGENT> sample(priorSampler());
//            endStats += sample.endState();
//        }
//        return BinomialDistribution(endStats);
//    }


//    double informationGain(const BinomialDistribution &analysis) {
//        ModelState<AGENT> realEndState = realTrajectory.endState();
//        BinomialDistribution prior = priorEndState(10000);
//        for(int i=0; i<realEndState.size(); ++i) {
//            std::cout << "Real occupancy = " << realEndState[i]
//            << " prior p = " << boost::math::pdf(prior.binomials[i], realEndState[i])
//            << " posterior p = " << boost::math::pdf(analysis.binomials[i], realEndState[i])
//            << " post / prior p = " << boost::math::pdf(analysis.binomials[i], realEndState[i])/ boost::math::pdf(prior.binomials[i], realEndState[i])
//            << std::endl;
////            << " prior = (" << prior.binomials[i].trials() << ", " << boost::math::mean(prior.binomials[i]) << ") "
////            << " posterior = (" << posterior.binomials[i].trials() << ", " << boost::math::mean(posterior.binomials[i]) << ") " << std::endl;
//        }
////        debug(
////                std::cout << "prior extendedLogProb = " << prior.logP(realEndState) << "Posterior extendedLogProb = " << posterior.logP(realEndState) << std::endl;
////                std::cout << "real end state = " << realEndState << std::endl;
////                std::cout << "prior = " << prior << std::endl;
////                std::cout << "posterior = " << posterior << std::endl;
////                );
//        return (analysis.logP(realEndState) - prior.logP(realEndState))/log(2.0);
//    }

//    AssimilationWindow(const Trajectory<AGENT> &realTrajectory,
//                       const PoissonState<AGENT> &priorStartState,
//                       std::vector<Observation<AGENT>> observations,
//                       int nSamples,
//                       int burnInSamples) :
//            realTrajectory(realTrajectory),
//            observations(observations) {
//        doAnalysis(priorStartState, nSamples, burnInSamples);
//    }


//    void doAnalysis(const PoissonState<AGENT> &priorStartState, int nSamples, int burnInSamples) {
//        PredPreyProblem<AGENT> abm(realTrajectory.nTimesteps(), observations, [&](const Trajectory<AGENT> &trajectory) {
//            return priorStartState.extendedLogProb(trajectory(0));
//        });
//        abm.advBasis();
//        abm.warmUp();
//        SimplexMCMC mcmc(abm, abm.logProbFunc());
//
//        Trajectory<AGENT> firstGuessTrajectory(realTrajectory.nTimesteps(), priorStartState.nextSample());
//        mcmc.setLPState(firstGuessTrajectory);
//        mcmc.findFeasibleStartPoint();
//        assert(mcmc.abmSanityChecks());
//
//        for(int burnIn=0; burnIn < burnInSamples; ++burnIn) { mcmc.nextSample(); }
//        for(int n=0; n<nSamples; ++n) {
//            mcmc.nextSample();
//            debug(if(nSamples<1000 || n%1000 == 1) {
//                assert(abm.isValidSolution(mcmc.X()));
//                std::cout << "Got nextSample " << mcmc.X() << std::endl;
//            });
//            const Trajectory<AGENT> &trajectory = reinterpret_cast<const Trajectory<AGENT> &>(mcmc.X());
//            analysis += trajectory.endState();
//        }
//        debug(std::cout
//                      << "infeasible/feasible: " << mcmc.infeasibleStatistics.nSamples *100.0/mcmc.feasibleStatistics.nSamples << "%" << std::endl
//                      << "Feasible nextSample statistics:" << std::endl << mcmc.feasibleStatistics << std::endl
//                      << "Infeasible nextSample statistics:" << std::endl << mcmc.infeasibleStatistics << std::endl;
//        );
//    }

};

// Gnuplot &operator<<(Gnuplot &gp, const AssimilationWindow<PredPreyAgent> &window);



#endif //GLPKTEST_ASSIMILATIONWINDOW_H
