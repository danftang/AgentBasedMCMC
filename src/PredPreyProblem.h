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
#include "agents/PredPreyAgent.h"
#include "ABMPrior.h"
#include "ABMPriorSampler.h"
#include "Likelihood.h"
#include "BernoulliStartState.h"
#include "PoissonStartState.h"
#include "TableauNormMinimiser.h"
#include "PredPreyTrajectory.h"
#include "Basis.h"

template<int GRIDSIZE,int NTIMESTEPS>
class PredPreyProblem {
public:
    typedef PredPreyAgent<GRIDSIZE> agent_type;
    typedef PredPreyTrajectory<GRIDSIZE,NTIMESTEPS> trajectory_type;
    typedef typename trajectory_type::value_type value_type;

    double                              pObserveIfPresent;
    double                              pPredator; // start state parameters
    double                              pPrey;
    trajectory_type                     realTrajectory;
    double                              kappa;
    std::vector<std::pair<State<agent_type>,int>> observations;
    Basis<trajectory_type>              basisObj;

//    PredPreyProblem(const std::string &filename): PredPreyProblem() {
//        std::ifstream probFile(filename);
//        if(!probFile.good()) throw("Can't open problem file " + filename);
//        boost::archive::binary_iarchive probArchive(probFile);
//        PredPreyProblem<GRIDSIZE> problem;
//        probArchive >> *this;
//    }
//
//
    PredPreyProblem()=default;


    PredPreyProblem(double pPredator, double pPrey, double pMakeObservation, double pObserveIfPresent, double kappa):
            pObserveIfPresent(pObserveIfPresent),
            pPredator(pPredator),
            pPrey(pPrey),
            realTrajectory(startState().priorSampler()()),
            kappa(kappa),
            observations(Likelihood<trajectory_type>::generateObservations(realTrajectory, pMakeObservation, pObserveIfPresent)),
            basisObj(posterior())
    {
    }


    auto startState() {
        return PoissonStartState<trajectory_type>([pPredator=pPredator, pPrey=pPrey](agent_type agent) {
            return agent.type() == PredPreyAgentBase::PREDATOR?pPredator:pPrey;
        });
//        return BernoulliStartState<trajectory_type>([pPredator=pPredator, pPrey=pPrey](agent_type agent) {
//            return agent.type() == PredPreyAgentBase::PREDATOR?pPredator:pPrey;
//        });
    };


    const Basis<trajectory_type> &basis() { return basisObj; }

    const std::vector<SparseVec<value_type>> &basisVectors() { return basisObj.basisVectors; }

    trajectory_type randomInitialSolution() { return basis().basisToDomain(randomBasisCoord()); }

    std::vector<value_type> randomBasisCoord() {
        std::vector<value_type> coord(basisObj.basisVectors.size());
        for(int i=0; i<coord.size(); ++i) {
            coord[i] = Random::nextBool(0.01);
        }
        return coord;
    }



    ABMPrior<trajectory_type> prior() {
        return ABMPrior<trajectory_type>(startState());
    };

    Likelihood<trajectory_type> likelihood() {
        return Likelihood<trajectory_type>(observations, pObserveIfPresent);
    }

    ConstrainedFactorisedDistribution<trajectory_type> posterior() {
        return prior() * likelihood();
    }
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
    friend std::ostream &operator <<(std::ostream &out, const PredPreyProblem<GRIDSIZE,NTIMESTEPS> &predPreyProblem) {
        out << "pObserveIfPresent: " << predPreyProblem.pObserveIfPresent << std::endl;
        out << "pPredator: " << predPreyProblem.pPredator << std::endl;
        out << "pPrey: " << predPreyProblem.pPrey << std::endl;
        out << "kappa: " << predPreyProblem.kappa << std::endl;
        out << "Real trajectory: " << predPreyProblem.realTrajectory << std::endl;
        out << "Observations: " << predPreyProblem.observations << std::endl;
        out << "Basis (nVectors x domainSize): " << predPreyProblem.basisObj.basisVectors.size() << " x " << predPreyProblem.basisObj.origin.size() << std::endl;
        out << "(Gridsize x Timesteps): " << GRIDSIZE << " x " << NTIMESTEPS << std::endl;
        return out;
    }
//
private:
    friend class boost::serialization::access;

//    template <typename Archive>
//    void load(Archive &ar, const unsigned int version) {
//        ar & pPredator & pPrey & kappa & realTrajectory & observations & tableau;
//    }
//
//    template <typename Archive>
//    void save(Archive &ar, const unsigned int version) const {
//        ar & pPredator & pPrey & kappa & realTrajectory & observations & tableau;
//    }

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & pObserveIfPresent & pPredator & pPrey & kappa & realTrajectory & observations & basisObj;
    }


//    BOOST_SERIALIZATION_SPLIT_MEMBER();

};


#endif //GLPKTEST_PREDPREYPROBLEM_H
