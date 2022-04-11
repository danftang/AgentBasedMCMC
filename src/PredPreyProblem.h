//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_PREDPREYPROBLEM_H
#define GLPKTEST_PREDPREYPROBLEM_H

#include <vector>
#include <boost/serialization/vector.hpp>
#include "Event.h"
#include "StlStream.h"
#include "StateTrajectory.h"
#include "AgentStateObservation.h"
#include "agents/PredPreyAgent.h"
#include "Prior.h"
#include "Likelihood.h"
#include "BernoulliStartState.h"
#include "TableauNormMinimiser.h"

template<int GRIDSIZE>
class PredPreyProblem {
public:
    double                              pPredator;
    double                              pPrey;
    double                              kappa;
//    double                              alpha;
    Prior<PredPreyAgent<GRIDSIZE>>      prior;
    Trajectory<PredPreyAgent<GRIDSIZE>> realTrajectory;
    Likelihood<PredPreyAgent<GRIDSIZE>> likelihood;
    WeightedFactoredConvexDistribution<ABM::occupation_type> posterior;
    TableauNormMinimiser<ABM::occupation_type> tableau;

    PredPreyProblem(): realTrajectory(0) {}

    PredPreyProblem(int nTimesteps, double pPredator, double pPrey, double pMakeObservation, double pObserveIfPresent, double kappa):
    pPredator(pPredator),
    pPrey(pPrey),
    kappa(kappa),
//    alpha(alpha),
    prior(nTimesteps, startStatePrior()),
    realTrajectory(prior.nextSample()),
    likelihood(realTrajectory, pMakeObservation, pObserveIfPresent),
    posterior(likelihood * prior),
    tableau(posterior.constraints) {
        tableau.findMinimalBasis();
    }


    auto startStatePrior() const {
//        return PoissonStartState<PredPreyAgent<GRIDSIZE>>([pPredator = pPredator, pPrey = pPrey](PredPreyAgent<GRIDSIZE> agent) {
//            return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
//        });
        return BernoulliStartState<PredPreyAgent<GRIDSIZE>>([pPredator = pPredator, pPrey = pPrey](PredPreyAgent<GRIDSIZE> agent) {
            return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
        });
    }


    int nTimesteps() const { return realTrajectory.nTimesteps(); }


    friend std::ostream &operator <<(std::ostream &out, const PredPreyProblem &predPreyProblem) {
        out << "pPredator: " << predPreyProblem.pPredator << std::endl;
        out << "pPrey: " << predPreyProblem.pPrey << std::endl;
        out << "kappa: " << predPreyProblem.kappa << std::endl;
        out << "(Gridsize x Timesteps): " << GRIDSIZE << " x " << predPreyProblem.nTimesteps() << std::endl;
        out << "Real trajectory: " << predPreyProblem.realTrajectory << std::endl;
        out << "Observations: " << predPreyProblem.likelihood.observations << std::endl;
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void load(Archive &ar, const unsigned int version) {
        std::vector<AgentStateObservation<PredPreyAgent<GRIDSIZE>>> observations;
        ar & pPredator & pPrey & kappa & realTrajectory & observations & tableau;
        prior = Prior<PredPreyAgent<GRIDSIZE>>(realTrajectory.nTimesteps(), startStatePrior());
        likelihood = Likelihood<PredPreyAgent<GRIDSIZE>>(observations);
        posterior = likelihood * prior;
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        ar & pPredator & pPrey & kappa & realTrajectory & likelihood.observations & tableau;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();


//    template <typename Archive>
//    void serialize(Archive &ar, const unsigned int version) {
//        ar & pPredator & pPrey & realTrajectory & likelihood.observations;
//    }

};


#endif //GLPKTEST_PREDPREYPROBLEM_H
