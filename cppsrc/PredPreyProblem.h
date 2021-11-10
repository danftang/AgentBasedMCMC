//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_PREDPREYPROBLEM_H
#define GLPKTEST_PREDPREYPROBLEM_H

#include <vector>
#include <boost/serialization/vector.hpp>
#include "glpkpp.h"
#include "Event.h"
#include "StlStream.h"
#include "StateTrajectory.h"
#include "AgentStateObservation.h"
#include "BernoulliModelState.h"
#include "agents/PredPreyAgent.h"

template<int GRIDSIZE>
class PredPreyProblem {
public:
    double pPredator;
    double pPrey;
    Trajectory<PredPreyAgent<GRIDSIZE>> realTrajectory;
    std::vector<AgentStateObservation<PredPreyAgent<GRIDSIZE>>> observations;

    PredPreyProblem(): realTrajectory(0) {}

    PredPreyProblem(int nTimesteps, double pPredator, double pPrey, double pMakeObservation, double pObserveIfPresent):
    pPredator(pPredator),
    pPrey(pPrey),
    realTrajectory(nTimesteps, startStatePrior().sampler()),
    observations(AgentStateObservation<PredPreyAgent<GRIDSIZE>>::generateObservations(realTrajectory, pMakeObservation, pObserveIfPresent)) {
    }


    BernoulliModelState<PredPreyAgent<GRIDSIZE>> startStatePrior() const {
        return BernoulliModelState<PredPreyAgent<GRIDSIZE>>([pPredator = pPredator, pPrey = pPrey](PredPreyAgent<GRIDSIZE> agent) {
            return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
        });
    }

    ConvexPMF<Trajectory<PredPreyAgent<GRIDSIZE>>> prior() const {
        return ConvexPMF<Trajectory<PredPreyAgent<GRIDSIZE>>>(nTimesteps(), startStatePrior().PMF());
    }

    std::function<Trajectory<PredPreyAgent<GRIDSIZE>>()> priorSampler() const {
        return Trajectory<PredPreyAgent<GRIDSIZE>>::priorSampler(nTimesteps(), startStatePrior().sampler());
    }

    ConvexPMF<Trajectory<PredPreyAgent<GRIDSIZE>>> likelihood() const {
        return ConvexPMF<Trajectory<PredPreyAgent<GRIDSIZE>>>::likelihood(nTimesteps(), observations);
    }

    ConvexPMF<Trajectory<PredPreyAgent<GRIDSIZE>>> posterior() const {
        return likelihood() * prior();
    }

    int nTimesteps() const { return realTrajectory.nTimesteps(); }

    friend std::ostream &operator <<(std::ostream &out, const PredPreyProblem &predPreyProblem) {
        out << "pPredator: " << predPreyProblem.pPredator << std::endl;
        out << "pPrey: " << predPreyProblem.pPrey << std::endl;
        out << "(Gridsize x Timesteps): " << GRIDSIZE << " x " << predPreyProblem.nTimesteps() << std::endl;
        out << "Real trajectory: " << predPreyProblem.realTrajectory << std::endl;
        out << "Observations: " << predPreyProblem.observations << std::endl;
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & pPredator & pPrey & realTrajectory & observations;
    }

};


#endif //GLPKTEST_PREDPREYPROBLEM_H
