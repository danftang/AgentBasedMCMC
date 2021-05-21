//
// Created by daniel on 18/05/2021.
//

#ifndef GLPKTEST_MODELSTATE_H
#define GLPKTEST_MODELSTATE_H

#include <map>
#include "Random.h"

template<typename AGENT>
class ModelState: public std::map<AGENT,double> {
public:
    using std::map<AGENT,double>::operator [];

    double operator[](const AGENT &agent) const {
        auto iter = this->find(agent);
        return iter==this->end()?0.0:iter->second;
    }


    ModelState &operator +=(const std::vector<AGENT> &stateVector) {
        for(const AGENT &agent: stateVector) {
            (*this)[agent] += 1.0;
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

};


#endif //GLPKTEST_MODELSTATE_H
