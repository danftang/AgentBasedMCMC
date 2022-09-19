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
                this->addFactor(
                        SparseFunction<std::pair<double,bool>, const Trajectory<AGENT> &>(
                                [state](const Trajectory<AGENT> &trajectory) {
                                    return widenedAgentMultinomial(state, trajectory);
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
        static double expMinusKappa = exp(-ABM::kappa);
        ABM::occupation_type l1StateOccupation(0);
        double newP(0.0);
        bool exactValue = true;
        std::vector<double> actPMF = state.agent.timestep(trajectory.temporaryPartialModelState(state.time, state.agent.neighbours()));
        for (int act = 0; act < AGENT::actDomainSize(); ++act) {
            int actOccupation = trajectory[Event<AGENT>(state.time, state.agent, act)];
            if (actOccupation != 0) {
                double pAct = actPMF[act];
                if (actOccupation < 0) { // negative occupation widening
//                        std::cout << "Widening due to negative occupation" << std::endl;
                    pAct = expMinusKappa;
                    actOccupation = -actOccupation;
                    exactValue = false;
                } else if (actPMF[act] == 0.0) {
//                        std::cout << "Widening due to impossible act" << std::endl;
                    pAct = expMinusKappa; // impossible act widening
                    exactValue = false;
                }
                newP += actOccupation*log(pAct) - lgamma(actOccupation + 1);
                l1StateOccupation += actOccupation; // sum of absolute values
            }
        }
        newP += lgamma(l1StateOccupation + 1); // Phi factorial
        return std::pair(newP,exactValue);
    }


    std::function<const Trajectory<AGENT> &()> sampler() const {
        return [sample = Trajectory<AGENT>(nTimesteps), startStateSampler = startStateSampler]() mutable -> const Trajectory<AGENT> & {
            ModelState<AGENT> t0State = startStateSampler();
            ModelState<AGENT> t1State;
            for (int t = 0; t < sample.nTimesteps(); ++t) {
                for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                    AGENT agent(agentId);
                    int nAgents = t0State[agentId];
                    for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                        sample[Event<AGENT>(t, agent, actId)] = 0;
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
        };
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
