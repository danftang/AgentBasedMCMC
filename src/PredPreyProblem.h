//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_PREDPREYPROBLEM_H
#define GLPKTEST_PREDPREYPROBLEM_H

#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "Event.h"
#include "StlStream.h"
#include "StateTrajectory.h"
#include "NoisyAgentStateObservation.h"
#include "agents/PredPreyAgent.h"
#include "Prior.h"
#include "Likelihood.h"
#include "BernoulliStartState.h"
#include "PoissonStartState.h"
#include "TableauNormMinimiser.h"

template<int GRIDSIZE>
class PredPreyProblem {
public:
    double                              pPredator;
    double                              pPrey;
    double                              kappa;
    Prior<PredPreyAgent<GRIDSIZE>>      prior;
    Trajectory<PredPreyAgent<GRIDSIZE>> realTrajectory;
    Likelihood<PredPreyAgent<GRIDSIZE>> likelihood;
    WeightedFactoredConvexDistribution<ABM::occupation_type> posterior;
    TableauNormMinimiser<ABM::occupation_type> tableau;

    PredPreyProblem(): realTrajectory(0) {}


    PredPreyProblem(const std::string &filename): PredPreyProblem() {
        std::ifstream probFile(filename);
        if(!probFile.good()) throw("Can't open problem file " + filename);
        boost::archive::binary_iarchive probArchive(probFile);
        PredPreyProblem<GRIDSIZE> problem;
        probArchive >> *this;
    }


    PredPreyProblem(int nTimesteps, double pPredator, double pPrey, double pMakeObservation, double pObserveIfPresent, double kappa):
    pPredator(pPredator),
    pPrey(pPrey),
    kappa(kappa),
    prior(nTimesteps, startStatePrior()),
    realTrajectory(prior.nextSample(true)),
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


    void save(const std::string &filename) {

    }


    friend std::ostream &operator <<(std::ostream &out, const PredPreyProblem &predPreyProblem) {
        out << "Real trajectory: " << predPreyProblem.realTrajectory << std::endl;
        out << "Observations: " << predPreyProblem.likelihood.noisyObservations << std::endl;
        out << "pPredator: " << predPreyProblem.pPredator << std::endl;
        out << "pPrey: " << predPreyProblem.pPrey << std::endl;
        out << "kappa: " << predPreyProblem.kappa << std::endl;
        out << "(Gridsize x Timesteps): " << GRIDSIZE << " x " << predPreyProblem.nTimesteps() << std::endl;
        out << "Basis mean L0 norm: " << predPreyProblem.tableau.meanColumnL0Norm() << std::endl;
        out << "Basis dimension (rows x cols): " << predPreyProblem.tableau.rows.size() << " x " << predPreyProblem.tableau.cols.size() << std::endl;
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void load(Archive &ar, const unsigned int version) {
        std::vector<NoisyAgentStateObservation<PredPreyAgent<GRIDSIZE>>> observations;
        ar & pPredator & pPrey & kappa & realTrajectory & observations & tableau;
        prior = Prior<PredPreyAgent<GRIDSIZE>>(realTrajectory.nTimesteps(), startStatePrior());
        likelihood = Likelihood<PredPreyAgent<GRIDSIZE>>(observations);
        posterior = likelihood * prior;
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        ar & pPredator & pPrey & kappa & realTrajectory & likelihood.noisyObservations & tableau;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();

};


#endif //GLPKTEST_PREDPREYPROBLEM_H
