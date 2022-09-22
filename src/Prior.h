// Represents the prior distribution over an Agent based model
//
// Created by daniel on 01/04/2022.
//

#ifndef ABMCMC_PRIOR_H
#define ABMCMC_PRIOR_H

#include "ABM.h"
#include "StartStateDistribution.h"
#include "ConstrainedFactorisedDistribution.h"
#include "PoissonStartState.h"
#include "TrajectoryDependencies.h"

// DOMAIN should be a domain over ABM trajectories that inplements
// index operators over Events and States
template<typename DOMAIN, typename AGENT  = typename DOMAIN::agent_type>
class Prior: public ConstrainedFactorisedDistribution<DOMAIN> {
public:
    std::function<const ModelState<AGENT> &()> startStateSampler;
    int           nTimesteps;

//    Prior(): nTimesteps(0) { }

    template<class STARTSTATE>
    Prior(int NTimesteps, const STARTSTATE &startState):
        nTimesteps(NTimesteps),
        startStateSampler(startState.modelStateSampler)
    {
        this->addConstraints(DOMAIN::constraints(nTimesteps));
        addMultinomials();
        (*this) *= startState;
    }


//    void addContinuityConstraints() {
//        for(int time = 1; time < nTimesteps; ++time) {
//            for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
//                SparseVec<ABM::coefficient_type> coefficients;
//                // outgoing edges
//                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
//                    coefficients.insert(DOMAIN::indexOf(Event<AGENT>(time, agentState, act)), 1.0);
//                }
//                // incoming edges
//                for (const Event<AGENT> &inEdge: State<AGENT>::incomingEventsByState[agentState]) {
//                    coefficients.insert(DOMAIN::indexOf(Event<AGENT>(time-1,inEdge.agent(),inEdge.act())), -1.0);
//                }
//                EqualityConstraint<ABM::occupation_type> constraint(coefficients, 0);
//                this->addConstraint(constraint);
//            }
//        }
//    }


    void addMultinomials() {
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                State<AGENT> state(time, AGENT(agentId));
//                std::vector<int> stateDependencies = DOMAIN::coefficients(state).indices; // this state occupation
////                std::set<int> actDependencies(stateDependencies.begin(), stateDependencies.end());
//                std::set<int> actDependencies;
//                for(int actId = 0; actId < AGENT::actDomainSize(); actId++) {
//                    actDependencies.insert(DOMAIN::indexOf(Event<AGENT>(state.time, state.agent, actId)));
//                }
//                for(const AGENT &neighbour :state.agent.neighbours()) {
//                    for(int neighbourDependency: DOMAIN::coefficients(State<AGENT>(time, neighbour)).indices) {
//                        actDependencies.insert(neighbourDependency);
//                    }
//                }
//                this->addFactor(
//                        SparseFunction<std::pair<double,bool>, const DOMAIN &>(
//                                [state](const DOMAIN &trajectory) {
//                                    return widenedAgentMultinomial(state, trajectory, false);
//                                },
//                                actDependencies
//                            )
//                        );

                for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                    Event<typename DOMAIN::agent_type> event(time, agentId, actId);
                    TrajectoryDependencies dependencies = DOMAIN::agent_type::eventProbDependencies(event);
                    dependencies.events.push_back(event);
                    this->addFactor(
                            [event](const DOMAIN &trajectory) {
                                return widenedEventFactor(event, trajectory);
                            },
                            dependencies.template indexDependencies<DOMAIN>()
                    );
                }

                this->addFactor(
                        [state](const DOMAIN &trajectory) {
                            return std::pair(lgamma(std::max(trajectory[state], 0) + 1), true);
                        },
                        DOMAIN::coefficients(state).indices
                );
            }
        }
    }


    // Turns an agent timestep, pi(a, \Psi), unction into a function
    // \prod_a \pi(a, \Psi, \psi)^T^t_{\psi a}/T^t_{\psi a}!
    // widened to decay (roughly) exponentially over negative occupations and zero probability actions
    static std::pair<double, bool> widenedAgentMultinomial(const State<AGENT> &state, const DOMAIN &trajectory, const bool includePhiFactorial) {
        ABM::occupation_type stateOccupation(0);
        double logP(0.0);
        bool exactValue = true;
        // std::vector<double> actPMF = state.agent.timestep(
//        const ModelState<AGENT> &neighbours = neighbourModelState(trajectory, state);
        for (int act = 0; act < AGENT::actDomainSize(); ++act) {
            Event<AGENT> event(state.time, state.agent, act);
            int actOccupation = trajectory[event];
            if (actOccupation != 0) {
                double logpAct;
                if (actOccupation < 0) { // negative occupation widening
                    logpAct = -ABM::kappa;
                    actOccupation = -actOccupation;
                    exactValue = false;
                } else {
                    logpAct = AGENT::logEventProb(event, trajectory);
                    if (logpAct == -INFINITY) {
                        logpAct = -ABM::kappa; // impossible act widening
                        exactValue = false;
                    }
                }
                logP += actOccupation * logpAct - lgamma(actOccupation + 1);
                stateOccupation += actOccupation;
            }
        }
        if(includePhiFactorial) logP += lgamma(std::max(stateOccupation, 0) + 1); // Phi factorial
        return std::pair(logP, exactValue);
    }

    static std::pair<double, bool> widenedEventFactor(const Event<AGENT> &event, const DOMAIN &trajectory) {
        double logP = 0.0;
        double logpAct;
        bool exactValue = true;
        int actOccupation = trajectory[event];
        if (actOccupation != 0) {
            if (actOccupation < 0) { // negative occupation widening
                logpAct = -ABM::kappa;
                actOccupation = -actOccupation;
                exactValue = false;
            } else {
                logpAct = AGENT::logEventProb(event, trajectory);
                if (logpAct == -INFINITY) {
                    logpAct = -ABM::kappa; // impossible act widening
                    exactValue = false;
                }
            }
            logP = actOccupation * logpAct - lgamma(actOccupation + 1);
        }
        return std::pair(logP, exactValue);
    }



    // fills only the neighbours of an agent-state, leaving the other elements undefined
    static const ModelState<AGENT> &neighbourModelState(const DOMAIN &trajectory, const State<AGENT> &agentState) {
        static thread_local ModelState<AGENT> temporaryModelState;
        for(AGENT neighbour: agentState.agent.eventProbDependencies())
            temporaryModelState[neighbour] = trajectory[State<AGENT>(agentState.time, neighbour)];
        return temporaryModelState;
    }



    std::function<const DOMAIN &()> sampler() const {
        return [sample = DOMAIN(nTimesteps), startStateSampler = startStateSampler]() mutable -> const DOMAIN & {
            ModelState<AGENT> t0State = startStateSampler();
            ModelState<AGENT> t1State;
            for (int t = 0; t < sample.nTimesteps(); ++t) {
                for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                    AGENT agent(agentId);
                    State<AGENT> state(t, agent);
                    int nAgents = t0State[agentId];
                    for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                        sample[Event<AGENT>(t, agent, actId)] = 0;
                    }
                    std::vector<double> actPMF = agentActPMF(state, t0State);
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
        };
    }

    static std::vector<double> agentActPMF(const State<AGENT> &state, const ModelState<AGENT> &modelState) {
        std::vector<double> pmf;
        pmf.reserve(AGENT::actDomainSize());
        for(int actId=0; actId < AGENT::actDomainSize(); ++actId) {
            Event<AGENT> event(state.time, state.agent, actId);
            pmf.push_back(exp(AGENT::logEventProb(event, modelState)));
        }
        return pmf;
    }

//    Trajectory<AGENT> nextSample() const {
//        static thread_local Trajectory<AGENT> sample(nTimesteps);
//        ModelState<AGENT> t0State = startState.sampler().nextSample();
//        ModelState<AGENT> t1State;
//        for (int t = 0; t < nTimesteps; ++t) {
//            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
//                AGENT agent(agentId);
//                int nAgents = t0State[agentId];
//                for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
//                    sample[Event<AGENT>(t, agent, actId)] = 0.0;
//                }
//                std::vector<double> actPMF = agent.timestep(t0State);
//                for (int a = 0; a < nAgents; ++a) {
//                    int nextAct = Random::nextIntFromDiscrete(actPMF);
//                    sample[Event<AGENT>(t, agent, nextAct)] += 1;
//                    t1State += agent.consequences(nextAct);
//                }
//            }
//            t0State.setToZero();
//            t0State.swap(t1State);
//        }
//        return sample;
//    }

//    friend std::ostream &operator <<(std::ostream &out, const Prior &prior) {
//        out << "Prior timesteps = " << prior.nTimesteps << "\nPrior ModelState =\n" << prior.startState << std::endl;
//        out << "Prior Trajectory constraints =\n" << prior.constraints;
//        return out;
//    }

};


#endif //ABMCMC_PRIOR_H
