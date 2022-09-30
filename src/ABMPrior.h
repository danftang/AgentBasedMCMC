// Represents the prior distribution over an Agent based model
//
// Created by daniel on 01/04/2022.
//

#ifndef ABMCMC_ABMPRIOR_H
#define ABMCMC_ABMPRIOR_H

#include "ABM.h"
#include "ConstrainedFactorisedDistribution.h"

// DOMAIN should be a domain over ABM trajectories that inplements
// index operators over Events and States
template<typename DOMAIN, typename AGENT  = typename DOMAIN::agent_type>
class ABMPrior: public ConstrainedFactorisedDistribution<DOMAIN> {
public:
//    std::function<const ModelState<AGENT> &()> startStateSampler;

    template<class STARTSTATE>
    ABMPrior(const STARTSTATE &startState)
//        :startStateSampler(startState._modelStateSampler)
    {
        this->addConstraints(DOMAIN::constraints());
        addMultinomials();
        (*this) *= startState;
    }


//    void addContinuityConstraints() {
//        for(int time = 1; time < nTimesteps; ++time) {
//            for(int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
//                SparseVec<ABM::coefficient_type> coefficients;
//                // outgoing edges
//                for (int act = 0; act < AGENT::actDomainSize; ++act) {
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
        for(int time = 0; time < DOMAIN::nTimesteps; ++time) {
            for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
                State<AGENT> state(time, AGENT(agentId));
//                std::vector<int> stateDependencies = DOMAIN::coefficients(state).indices; // this state occupation
////                std::set<int> actDependencies(stateDependencies.begin(), stateDependencies.end());
//                std::set<int> actDependencies;
//                for(int actId = 0; actId < AGENT::actDomainSize; actId++) {
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

                for (int actId = 0; actId < AGENT::actDomainSize; ++actId) {
                    Event<typename DOMAIN::agent_type> event(time, agentId, actId);
                    std::vector<int> dependencies = DOMAIN::agent_type::template eventProbDependencies<DOMAIN>(event);
                    dependencies.push_back(DOMAIN::indexOf(event));
                    this->addFactor(
                            [event](const DOMAIN &trajectory) {
                                return widenedEventFactor(event, trajectory);
                            },
                            dependencies
                    );
                }

                this->addFactor( // factorial of state occupation
                        [state](const DOMAIN &trajectory) {
                            auto occupation = trajectory[state];
                            return (occupation < 0) ? std::pair(ABM::kappa*occupation,false) : std::pair(lgamma(occupation + 1), true); // TODO: Test!!
//                            return std::pair(lgamma(std::abs(trajectory[state]) + 1), true); // TODO: test
//                            return std::pair(lgamma(std::max(trajectory[state], 0) + 1), true);
                        },
                        {DOMAIN::indexOf(state)}
                );
            }
        }
    }


    // Turns an agent timestep, pi(a, \Psi), unction into a function
    // \prod_a \pi(a, \Psi, \psi)^T^t_{\psi a}/T^t_{\psi a}!
    // widened to decay (roughly) exponentially over negative occupations and zero probability actions
//    static std::pair<double, bool> widenedAgentMultinomial(const State<AGENT> &state, const DOMAIN &trajectory, const bool includePhiFactorial) {
//        ABM::occupation_type stateOccupation(0);
//        double logP(0.0);
//        bool exactValue = true;
//        // std::vector<double> actPMF = state.agent.timestep(
////        const ModelState<AGENT> &neighbours = neighbourModelState(trajectory, state);
//        for (int act = 0; act < AGENT::actDomainSize; ++act) {
//            Event<AGENT> event(state.time, state.agent, act);
//            int actOccupation = trajectory[event];
//            if (actOccupation != 0) {
//                double logpAct;
//                if (actOccupation < 0) { // negative occupation widening
//                    logpAct = -ABM::kappa;
//                    actOccupation = -actOccupation;
//                    exactValue = false;
//                } else {
//                    logpAct = AGENT::logEventProb(event, trajectory);
//                    if (logpAct == -INFINITY) {
//                        logpAct = -ABM::kappa; // impossible act widening
//                        exactValue = false;
//                    }
//                }
//                logP += actOccupation * logpAct - lgamma(actOccupation + 1);
//                stateOccupation += actOccupation;
//            }
//        }
//        if(includePhiFactorial) logP += lgamma(std::max(stateOccupation, 0) + 1); // Phi factorial
//        return std::pair(logP, exactValue);
//    }

    static std::pair<double, bool> widenedEventFactor(const Event<AGENT> &event, const DOMAIN &trajectory) {
        int actOccupation = trajectory[event];
        if(actOccupation == 0) return std::pair(0.0, true);
        if (actOccupation < 0) { // negative occupation widening
//            return std::pair(-actOccupation * log(1.0/AGENT::actDomainSize) - lgamma(1-actOccupation) + ABM::kappa * actOccupation, false); // pair production and decay
//            return std::pair(ABM::kappa * actOccupation - lgamma(1-actOccupation), false);
            return std::pair(ABM::kappa * actOccupation, false);
        }
        auto [logpAct, isFeasible] = AGENT::widenedLogEventProb(event, trajectory); // TODO: test!!!
//        if (actOccupation > 1) { // Fermionic occupation widening
//            return std::pair(logpAct + ABM::kappa * (1-actOccupation), false); // TODO: test!!!
//        }
        return std::pair(actOccupation * logpAct - lgamma(actOccupation + 1), isFeasible);

//        double logpAct = AGENT::logEventProb(event, trajectory);
//        if (logpAct == -INFINITY) { // impossible act widening
//            return std::pair(event.agent().logMarginalTimestep(event.act())* actOccupation - lgamma(actOccupation + 1) - ABM::kappa, false);
//        }
//        return std::pair(actOccupation * logpAct - lgamma(actOccupation + 1), true);
    }



    // fills only the neighbours of an agent-state, leaving the other elements undefined
    static const ModelState<AGENT> &neighbourModelState(const DOMAIN &trajectory, const State<AGENT> &agentState) {
        static thread_local ModelState<AGENT> temporaryModelState;
        for(AGENT neighbour: agentState.agent.eventProbDependencies())
            temporaryModelState[neighbour] = trajectory[State<AGENT>(agentState.time, neighbour)];
        return temporaryModelState;
    }



//    std::function<const DOMAIN &()> sampler() const {
//        return [sample = DOMAIN(), startStateSampler = startStateSampler]() mutable -> const DOMAIN & {
//            sample.setStartState(startStateSampler());
//            for (int t = 0; t < DOMAIN::nTimesteps; ++t) {
//                if(t>0) sample.recalculateDependentVariables(t);
//                for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
//                    AGENT agent(agentId);
//                    State<AGENT> state(t, agent);
//                    int nAgents = sample[state];
//                    for (int actId = 0; actId < AGENT::actDomainSize; ++actId) {
//                        sample[DOMAIN::indexOf(Event<AGENT>(t, agent, actId))] = 0; // modification must be via int index
//                    }
//                    std::vector<double> actPMF = agentActPMF(state, sample);
//                    for (int a = 0; a < nAgents; ++a) {
//                        int nextAct = Random::nextIntFromDiscrete(actPMF);
//                        sample[DOMAIN::indexOf(Event<AGENT>(t, agent, nextAct))] += 1;
//                    }
//                }
//            }
//            return sample;
//        };
//    }


//    Trajectory<AGENT> nextSample() const {
//        static thread_local Trajectory<AGENT> sample(nTimesteps);
//        ModelState<AGENT> t0State = startState.sampler().nextSample();
//        ModelState<AGENT> t1State;
//        for (int t = 0; t < nTimesteps; ++t) {
//            for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
//                AGENT agent(agentId);
//                int nAgents = t0State[agentId];
//                for (int actId = 0; actId < AGENT::actDomainSize; ++actId) {
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


#endif //ABMCMC_ABMPRIOR_H
