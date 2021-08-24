//
// Created by daniel on 18/05/2021.
//

#ifndef GLPKTEST_MODELSTATE_H
#define GLPKTEST_MODELSTATE_H

#include <map>
#include "Random.h"

// represents occupation numbers of agent states in a single timestep.
template<typename AGENT>
class ModelState: public std::vector<double> {
public:

    ModelState(): std::vector<double>(AGENT::domainSize(),0.0) { };

    ModelState(std::vector<double> &&X): std::vector<double>(X) {
        assert(size() == AGENT::domainSize());
    }

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
        for(int i=0; i<size(); ++i) (*this)[i] = 0.0;
    }

    void clear() {
        std::cerr << "ERROR: don't try to resize a modelstate!" << std::endl;
        assert(false);
    }

    void resize(size_type i) {
        std::cerr << "ERROR: don't try to resize a modelstate!" << std::endl;
        assert(false);
    }

    void resize(size_type i, double d) {
        std::cerr << "ERROR: don't try to resize a modelstate!" << std::endl;
        assert(false);
    }


    ModelState &operator +=(const std::vector<AGENT> &stateVector) {
        for(const AGENT &agent: stateVector) {
            (*this)[agent] += 1.0;
        }
        return *this;
    }


    // Addition for aggregation of states
    ModelState &operator +=(const ModelState<AGENT> &other) {
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            (*this)[agentId] += other[agentId];
        }
        return *this;
    }

    ModelState &operator *=(double multiplier) {
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            (*this)[agentId] *= multiplier;
        }
        return *this;
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
