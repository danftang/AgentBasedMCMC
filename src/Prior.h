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

template<class AGENT>
class Prior: public ConstrainedFactorisedDistribution<const Trajectory<AGENT> &, ABM::coefficient_type> {
public:
    PoissonStartState<AGENT>       startState;
    int                            nTimesteps;

    using ConstrainedFactorisedDistribution<const Trajectory<AGENT> &,ABM::coefficient_type>::addFactor;

    Prior(): nTimesteps(0) { }

    Prior(int NTimesteps, PoissonStartState<AGENT> StartState):
        nTimesteps(NTimesteps),
        startState(std::move(StartState))
    {
        addContinuityConstraints();
        addInteractionFactors();
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
                for(const AGENT &neighbour :state.agent.dependencies()) {
                    for(int neighbourAct: State<AGENT>(time, neighbour).forwardOccupationDependencies()) {
                        actDependencies.push_back(neighbourAct);
                    }
                }
                addFactor(
                        SparseWidenedFunction<double, Trajectory<AGENT>>(
                                [state](const Trajectory<AGENT> &trajectory) {
                                    return Prior<AGENT>::widenedAgentMultinomial(state, trajectory);
                                },
                                state.forwardOccupationDependencies()
                            ),
                            actDependencies
                        );
            }
        }
    }


    // Turns an agent timestep, pi(a, \Psi), unction into a function
    // \prod_a \pi(a, \Psi, \psi)^T^t_{\psi a}/T^t_{\psi a}!
    // widened to decay exponentially over negative occupations and zero probability actions
    static std::pair<double, bool> widenedAgentMultinomial(const State<AGENT> &state, const Trajectory<AGENT> &trajectory) {
        ABM::occupation_type stateOccupation = trajectory[state];
        double newLogP = (stateOccupation > 1)?lgamma(stateOccupation + 1.0):0.0; // log of Phi factorial
        bool exactValue = true;
        if(stateOccupation > 0) {
            std::vector<double> actPMF = state.agent.timestep(trajectory);
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
        // TOSO: implement this
    }

};


#endif //ABMCMC_PRIOR_H
