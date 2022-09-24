//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_PREDPREYPROBLEM_H
#define GLPKTEST_PREDPREYPROBLEM_H

#include <vector>
#include <fstream>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "Event.h"
#include "include/StlStream.h"
#include "StateTrajectory.h"
#include "agents/PredPreyAgent.h"
#include "Prior.h"
#include "Likelihood.h"
#include "BernoulliStartState.h"
#include "PoissonStartState.h"
#include "TableauNormMinimiser.h"

template<int GRIDSIZE>
class PredPreyProblem {
public:
//    typedef PredPreyAgent<GRIDSIZE> agent_type;
//    typedef ExtendedTrajectory2<GRIDSIZE,> trajectory_type;
//
//    double                              pObserveIfPresent;
//    std::vector<std::pair<State<agent_type>,int>> observations;
//
//    double                              pPredator; // start state parameters
//    double                              pPrey;
//
//    trajectory_type realTrajectory;
//
//    double                              kappa;
////    Prior<PoissonStartState<agent_type>>      prior;
////    Likelihood<agent_type> likelihood;
////    ConstrainedFactorisedDistribution<Trajectory<agent_type>> posterior;
//    TableauNormMinimiser<ABM::occupation_type> tableau;
//
//    PredPreyProblem(): realTrajectory(0) {}
//
////    PredPreyProblem(const std::string &filename): PredPreyProblem() {
////        std::ifstream probFile(filename);
////        if(!probFile.good()) throw("Can't open problem file " + filename);
////        boost::archive::binary_iarchive probArchive(probFile);
////        PredPreyProblem<GRIDSIZE> problem;
////        probArchive >> *this;
////    }
//
//
//    PredPreyProblem(int nTimesteps, double pPredator, double pPrey, double pMakeObservation, double pObserveIfPresent, double kappa):
//    pObserveIfPresent(pObserveIfPresent),
//    pPredator(pPredator),
//    pPrey(pPrey),
//    realTrajectory(prior().sampler()()),
//    kappa(kappa),
//    observations(Likelihood<trajectory_type>::generateObservaions(realTrajectory, pMakeObservation, pObserveIfPresent))
//    {
//
//    }
//
//    Prior<PoissonStartState<agent_type>> prior() {
//        return Prior(nTimesteps, startStatePrior());
//    };
//
//    Likelihood<trajectory_type> likelihood() {
//        return Likelihood<trajectory_type>(observations, pObserveIfPresent);
//    }
//
//
//    ConstrainedFactorisedDistribution<Trajectory<agent_type>> posterior() {
//        return prior() * likelihood();
//    }
//
//
//    auto startStatePrior() const {
//        return PoissonStartState<PredPreyAgent<GRIDSIZE>>([pPredator = pPredator, pPrey = pPrey](PredPreyAgent<GRIDSIZE> agent) {
//            return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
//        });
////        return BernoulliStartState<PredPreyAgent<GRIDSIZE>>([pPredator = pPredator, pPrey = pPrey](PredPreyAgent<GRIDSIZE> agent) {
////            return agent.type() == PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
////        });
//    }
//
//
//    int nTimesteps() const { return realTrajectory.nTimesteps(); }
//
//
//    friend std::ostream &operator <<(std::ostream &out, const PredPreyProblem &predPreyProblem) {
//        out << "Real trajectory: " << predPreyProblem.realTrajectory << std::endl;
//        out << "Observations: " << predPreyProblem.likelihood.noisyObservations << std::endl;
//        out << "pPredator: " << predPreyProblem.pPredator << std::endl;
//        out << "pPrey: " << predPreyProblem.pPrey << std::endl;
//        out << "kappa: " << predPreyProblem.kappa << std::endl;
//        out << "(Gridsize x Timesteps): " << GRIDSIZE << " x " << predPreyProblem.nTimesteps() << std::endl;
//        out << "Basis mean L0 norm: " << predPreyProblem.tableau.meanColumnL0Norm() << std::endl;
//        out << "Basis dimension (rows x cols): " << predPreyProblem.tableau.rows.size() << " x " << predPreyProblem.tableau.cols.size() << std::endl;
//        return out;
//    }
//
//private:
//    friend class boost::serialization::access;
//
//    template <typename Archive>
//    void load(Archive &ar, const unsigned int version) {
//        std::vector<NoisyAgentStateObservation<PredPreyAgent<GRIDSIZE>>> observations;
//        ar & pPredator & pPrey & kappa & realTrajectory & observations & tableau;
//    }
//
//    template <typename Archive>
//    void save(Archive &ar, const unsigned int version) const {
//        ar & pPredator & pPrey & kappa & realTrajectory & observations & tableau;
//    }
//
//    BOOST_SERIALIZATION_SPLIT_MEMBER();

};


#endif //GLPKTEST_PREDPREYPROBLEM_H
