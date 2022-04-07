//
// Created by daniel on 18/05/2021.
//

#ifndef GLPKTEST_MODELSTATE_H
#define GLPKTEST_MODELSTATE_H

#include <map>
#include <functional>
#include <valarray>
#include "Random.h"
#include "ABM.h"
#include "State.h"
// #include "Trajectory.h"

// represents occupation numbers of agent states in a single timestep.
template<typename AGENT>
class ModelState: public std::vector<ABM::occupation_type> {
public:

    ModelState(): std::vector<value_type>(AGENT::domainSize(),0) { };

    ModelState(std::vector<value_type> &&X): std::vector<value_type>(X) {
        assert(size() == AGENT::domainSize());
    }


    ModelState(const std::vector<ABM::occupation_type> &trajectory, int nTimesteps, int time) : ModelState() {
        assert(time>=0 && time<=nTimesteps);
        if(time < nTimesteps) {
            for (int stateId = 0; stateId < AGENT::domainSize(); ++stateId)
                (*this)[stateId] = State<AGENT>(time,stateId).forwardOccupation(trajectory);
        } else {
            for (int stateId = 0; stateId < AGENT::domainSize(); ++stateId)
                (*this)[stateId] = State<AGENT>(time,stateId).backwardOccupation(trajectory);
        }
    }

//    ModelState(const Trajectory<AGENT> &trajectory, int time): ModelState(trajectory, trajectory.nTimesteps(), time) {}


    // ensure zero initialisation
//    double &operator[](const AGENT &agent) {
//        auto iter = this->find(agent);
//        if(iter == this->end()) {
//            double &newEntry = std::map<AGENT,double>::operator[](agent);
//            newEntry = 0.0;
//            return newEntry;
//        }
//        return iter->second;
//    }
//
//    double operator[](const AGENT &agent) const {
//        auto iter = this->find(agent);
//        return iter==this->end()?0.0:iter->second;
//    }

    void setToZero() {
        for(int i=0; i<size(); ++i) (*this)[i] = 0;
    }

    void clear() {
        std::cerr << "ERROR: don't try to resize a modelstate!" << std::endl;
        assert(false);
    }

    void resize(size_t i) {
        std::cerr << "ERROR: don't try to resize a modelstate!" << std::endl;
        assert(false);
    }

    void resize(size_t i, double d) {
        std::cerr << "ERROR: don't try to resize a modelstate!" << std::endl;
        assert(false);
    }


    ModelState<AGENT> &operator +=(const std::vector<AGENT> &stateVector) {
        for(const AGENT &agent: stateVector) {
            (*this)[agent] += 1;
        }
        return *this;
    }


    // Addition for aggregation of states
    ModelState<AGENT> &operator +=(const ModelState<AGENT> &other) {
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            (*this)[agentId] += other[agentId];
        }
        return *this;
    }

    ModelState<AGENT> &operator *=(value_type multiplier) {
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            (*this)[agentId] *= multiplier;
        }
        return *this;
    }

    std::valarray<double> operator /(double denominator) const {
        std::valarray<double> probs(size());
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            probs[agentId] = (*this)[agentId] / denominator;
        }
        return probs;
    }

    std::valarray<double> operator *(double multiplier) const {
        std::valarray<double> probs(size());
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            probs[agentId] = (*this)[agentId] * multiplier;
        }
        return probs;
    }


    static ModelState<AGENT> randomPoissonState(const std::function<double (const AGENT &)> &pmf) {
        ModelState<AGENT> state;
        for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
            int occupation = Random::nextPoisson(pmf(agentId));
            if(occupation > 0) state[agentId] = occupation;
        }
        return state;
    }

//    friend std::ostream &operator <<(std::ostream &out, const ModelState<AGENT> &modelState) {
//        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
//            if(modelState[agentId] != 0.0) out << AGENT(agentId) << " -> " << modelState[agentId] << " ";
//        }
//        return  out;
//    }

};

#endif //GLPKTEST_MODELSTATE_H
