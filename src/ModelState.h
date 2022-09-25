// represents occupation numbers of agent states in a single timestep.
//
// Created by daniel on 18/05/2021.
//

#ifndef GLPKTEST_MODELSTATE_H
#define GLPKTEST_MODELSTATE_H

#include <map>
#include <functional>
#include <valarray>
#include "include/Random.h"
#include "ABM.h"
#include "State.h"

template<typename AGENT>
class ModelState: public std::vector<int> {
public:
    typedef AGENT agent_type;

    static const ModelState<AGENT> zero;

    ModelState(): std::vector<value_type>(AGENT::domainSize,0) { };

    ModelState(std::vector<value_type> &&X): std::vector<value_type>(X) {
        assert(size() == AGENT::domainSize);
    }

    template<class DOMAIN>
    ModelState(const DOMAIN &trajectory, int time);

    // so that Modelstate can also look like a trajectory
    const value_type &operator[](const State<AGENT> &state) const {
        return std::vector<value_type>::operator[]((int)state.agent);
    }

    const value_type &operator[](const AGENT &agent) const {
        return std::vector<value_type>::operator[]((int)agent);
    }

    value_type &operator[](const AGENT &agent) {
        return std::vector<value_type>::operator[]((int)agent);
    }

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
        for(int agentId=0; agentId < AGENT::domainSize; ++agentId) {
            (*this)[agentId] += other[agentId];
        }
        return *this;
    }

    ModelState<AGENT> &operator *=(value_type multiplier) {
        for(int agentId=0; agentId < AGENT::domainSize; ++agentId) {
            (*this)[agentId] *= multiplier;
        }
        return *this;
    }

    std::valarray<double> operator /(double denominator) const {
        std::valarray<double> probs(size());
        for(int agentId=0; agentId < AGENT::domainSize; ++agentId) {
            probs[agentId] = (*this)[agentId] / denominator;
        }
        return probs;
    }

    std::valarray<double> operator *(double multiplier) const {
        std::valarray<double> probs(size());
        for(int agentId=0; agentId < AGENT::domainSize; ++agentId) {
            probs[agentId] = (*this)[agentId] * multiplier;
        }
        return probs;
    }


//    static ModelState<AGENT> randomPoissonState(const std::function<double (const AGENT &)> &pmf) {
//        ModelState<AGENT> state;
//        for(int agentId=0; agentId<AGENT::domainSize; ++agentId) {
//            int occupation = Random::nextPoisson(pmf(agentId));
//            if(occupation > 0) state[agentId] = occupation;
//        }
//        return state;
//    }

};

template<class AGENT>
template<class DOMAIN>
ModelState<AGENT>::ModelState(const DOMAIN &trajectory, int time) : ModelState() {
    for (int stateId = 0; stateId < AGENT::domainSize; ++stateId)
        (*this)[stateId] = trajectory[State<AGENT>(time,stateId)];
}

template<typename AGENT> const ModelState<AGENT> ModelState<AGENT>::zero;


#endif //GLPKTEST_MODELSTATE_H
