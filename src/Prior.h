// Represents the prior distribution over an ABM
//
// Created by daniel on 01/04/2022.
//

#ifndef ABMCMC_PRIOR_H
#define ABMCMC_PRIOR_H

#include "ABM.h"
#include "StartStateDistribution.h"
#include "ConstrainedFactorisedDistribution.h"
#include "PoissonStartState.h"

template<typename STARTSTATE, typename AGENT = typename STARTSTATE::domain_type::agent_type>
class Prior: public ConstrainedFactorisedDistribution<Trajectory<AGENT>, ABM::coefficient_type> {
public:
    STARTSTATE       startState;
    int              nTimesteps;

    using ConstrainedFactorisedDistribution<Trajectory<AGENT>,ABM::coefficient_type>::addFactor;

    Prior(): nTimesteps(0) { }

    Prior(int NTimesteps, STARTSTATE StartState):
        nTimesteps(NTimesteps),
        startState(std::move(StartState))
    {
        addContinuityConstraints();
        addInteractionFactors();
        (*this) *= startState;
    }


    void addContinuityConstraints() {
        for(int time = 1; time < nTimesteps; ++time) {
            for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                SparseVec<ABM::coefficient_type> coefficients;
                // outgoing edges
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    coefficients.insert(Event<AGENT>(time, agentState, act).id, 1.0);
                }
                // incoming edges
                for (const Event<AGENT> &inEdge: State<AGENT>::incomingEventsByState[agentState]) {
                    coefficients.insert(Event<AGENT>(time-1,inEdge.agent(),inEdge.act()).id, -1.0);
                }
                EqualityConstraint<ABM::occupation_type> constraint(coefficients, 0);
                this->addConstraint(constraint);
            }
        }
    }

    void addInteractionFactors() {
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                State<AGENT> state(time, AGENT(agentId));
                std::vector<int> actDependencies = state.forwardOccupationDependencies();
                for(const AGENT &neighbour :state.agent.neighbours()) {
                    for(int neighbourAct: State<AGENT>(time, neighbour).forwardOccupationDependencies()) {
                        actDependencies.push_back(neighbourAct);
                    }
                }
                addFactor(
                        SparseFunction<std::pair<double,bool>, const Trajectory<AGENT> &>(
                                [state](const Trajectory<AGENT> &trajectory) {
                                    return Prior<STARTSTATE,AGENT>::widenedAgentMultinomial(state, trajectory);
                                },
                                actDependencies
                            )
                        );
            }
        }
    }


    // Turns an agent timestep, pi(a, \Psi), unction into a function
    // \prod_a \pi(a, \Psi, \psi)^T^t_{\psi a}/T^t_{\psi a}!
    // widened to decay (roughly) exponentially over negative occupations and zero probability actions
    static std::pair<double, bool> widenedAgentMultinomial(const State<AGENT> &state, const Trajectory<AGENT> &trajectory) {
        ABM::occupation_type stateOccupation = trajectory[state];
        double newLogP = (stateOccupation > 1)?lgamma(stateOccupation + 1.0):0.0; // log of Phi factorial
        bool exactValue = true;
        if(stateOccupation > 0) {
            std::vector<double> actPMF = state.agent.timestep(trajectory, state.time);
            for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                int actOccupation = trajectory[Event<AGENT>(state.time, state.agent, act)];
                double pAct = actPMF[act];
                if(actPMF[act] == 0.0) {
                    pAct = exp(-ABM::kappa); // impossible act widening
                    exactValue = false;
                }
                if(actOccupation < 0.0) { // negative occupation widening
                    pAct = exp(-ABM::kappa);
                    actOccupation = -actOccupation;
                    exactValue = false;
                }
                newLogP += actOccupation*log(pAct) - lgamma(actOccupation + 1.0); // Fermionic bounding
            }
        }
        return std::pair(newLogP,exactValue);
    }


    Trajectory<AGENT> nextSample() const {
        Trajectory<AGENT> sample(nTimesteps);
        ModelState<AGENT> t0State = startState.nextSample();
        ModelState<AGENT> t1State;
        for (int t = 0; t < nTimesteps; ++t) {
            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                AGENT agent(agentId);
                int nAgents = t0State[agentId];
                for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                    sample[Event<AGENT>(t, agent, actId)] = 0.0;
                }
                std::vector<double> actPMF = agent.timestep(t0State);
                for (int a = 0; a < nAgents; ++a) {
                    int nextAct = Random::nextIntFromDiscrete(actPMF);
                    sample[Event<AGENT>(t, agent, nextAct)] += 1;
                    t1State += agent.consequences(nextAct);
                }
            }
            t0State.setToZero();
            t0State.swap(t1State);
        }
        return sample;
    }

    friend std::ostream &operator <<(std::ostream &out, const Prior &prior) {
        out << "Prior timesteps = " << prior.nTimesteps << " Start state =\n" << prior.startState << std::endl;
        out << prior.constraints;
        return out;
    }

};


#endif //ABMCMC_PRIOR_H
