//
// Created by daniel on 26/09/22.
//

#ifndef ABMCMC_ABMPRIORSAMPLER_H
#define ABMCMC_ABMPRIORSAMPLER_H

#include <functional>
#include "ModelState.h"

template<class DOMAIN, class AGENT = typename DOMAIN::agent_type>
class ABMPriorSampler {
public:
    DOMAIN sample;
    std::function<const ModelState<AGENT> &()> startStateSampler;

    ABMPriorSampler(std::function<const ModelState<AGENT> &()> startStateSampler): startStateSampler(std::move(startStateSampler)) {}

    const DOMAIN &operator()() {
        sample.setStartState(startStateSampler());
        for (int t = 0; t < DOMAIN::nTimesteps; ++t) {
            if (t > 0) sample.recalculateDependentVariables(t);
            for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
                AGENT agent(agentId);
                State<AGENT> state(t, agent);
                int nAgents = sample[state];
                for (int actId = 0; actId < AGENT::actDomainSize; ++actId) {
                    sample[DOMAIN::indexOf(Event<AGENT>(t, agent, actId))] = 0; // modification must be via int index
                }
                std::vector<double> actPMF = agentActPMF(state, sample);
                for (int a = 0; a < nAgents; ++a) {
                    int nextAct = Random::nextIntFromDiscrete(actPMF);
                    sample[DOMAIN::indexOf(Event<AGENT>(t, agent, nextAct))] += 1;
                }
            }
        }
        return sample;
    }

    static std::vector<double> agentActPMF(const State<AGENT> &state, const DOMAIN &trajectory) {
        std::vector<double> pmf;
        pmf.reserve(AGENT::actDomainSize);
        for(int actId=0; actId < AGENT::actDomainSize; ++actId) {
            Event<AGENT> event(state.time, state.agent, actId);
            pmf.push_back(exp(AGENT::logEventProb(event, trajectory)));
        }
        return pmf;
    }

};


#endif //ABMCMC_ABMPRIORSAMPLER_H
