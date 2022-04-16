// Represents the probability distribution of a trajectory given the start state,
// where the start state is implicit in the trajectory.
//
// The distribution is expressed as a weighted, factored convex distribution consisting of
// three sets of factors
// 1) the act fermionic factors, of the form
// 0 <= X_i <= 1
// F_i(X_i) = E_{P(Phi|psi(X_i)=1)}[pi(psi(X_i), Phi, a(X_i))]
//
//
//
// 2) the continuity constraints
// 3) the interaction constraints
//
// ...and a weight given by
//
// \prod_psi Phi_psi! \prod_a (pi(psi, Phi, a) / pi(psi, PhiBar, a))^{T^t_{psi a}}
//
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_FORECAST_H
#define ABMCMC_FORECAST_H


#include "WeightedFactoredConvexDistribution.h"
#include "TrajectoryImportance.h"
#include "ABM.h"
#include "Event.h"
#include "StartStateDistribution.h"

template <class AGENT>
class Forecast: public WeightedFactoredConvexDistribution<ABM::occupation_type> {
public:
    int nTimesteps;

    explicit Forecast(int NTimesteps): nTimesteps(NTimesteps),
                                       WeightedFactoredConvexDistribution<ABM::occupation_type>( [NTimesteps]() {
        return std::unique_ptr<PerturbableFunction<ABM::occupation_type,double>>(new TrajectoryImportance<AGENT>(NTimesteps));
    }) {
        addActFermionicFactors();
        addContinuityConstraints();
        addInteractionConstraints();
    }


    void addContinuityConstraints() {
        for(int time = 1; time < nTimesteps; ++time) {
            for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                Constraint<ABM::occupation_type> constraint(0.0 ,0.0);
                // outgoing edges
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    constraint.coefficients.insert(Event<AGENT>(time, agentState, act).id, 1.0);
                }
                // incoming edges
                for (const Event<AGENT> &inEdge: State<AGENT>::incomingEventsByState[agentState]) {
                    constraint.coefficients.insert(Event<AGENT>(time-1,inEdge.agent(),inEdge.act()).id, -1.0);
                }
                this->addFactor(constraint);
            }
        }
    }


    void addInteractionConstraints() {
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                AGENT agent(agentState);
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    for(const Constraint<ABM::occupation_type> &actConstraint : agent.constraints(time, act)) {
                        push_xImpliesY(Event<AGENT>(time,agent,act).id, actConstraint);
                    }
                }
            }
        }
    }


    void addActFermionicFactors() {
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    Event<AGENT> event = Event<AGENT>(time,agentState,act);
                    double logMarginalPi = log(event.agent().marginalTimestep(event.act()));
                    this->addFactor(
                            0 <= 1*event <= 1,
                            [logMarginalPi](ABM::occupation_type occupation) {
                                return occupation*logMarginalPi;
                            }
                    );
                }
            }
        }
    }

    Trajectory<AGENT> nextSample(const StartStateDistribution<AGENT> &startStateDist, bool isFermionic) const {
        return nextSample(startStateDist, nTimesteps, isFermionic);
    }

    // assumes act-Fermionicity
    static Trajectory<AGENT> nextSample(const StartStateDistribution<AGENT> &startStateDist, int nTimesteps, bool isFermionic) {
        Trajectory<AGENT> sample(nTimesteps);
        bool isValid;
        int nAttempts = 0;
        do {
            ModelState<AGENT> t0State = startStateDist.nextSample();
            ModelState<AGENT> t1State;
            isValid = true;
            for (int t = 0; t < nTimesteps; ++t) {
                for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                    AGENT agent(agentId);
                    int nAgents = t0State[agentId];
                    for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                        sample[Event<AGENT>(t, agent, actId)] = 0.0;
                    }
                    // now choose acts for each of nAgents from act Fermionic distribution

                    std::vector<double> actPMF = agent.timestep(t0State);
//                    std::cout << "Got act distribution " << actPMF << std::endl;

                    if(isFermionic) {
                        std::vector<bool> chosenActs(actPMF.size(), false);
                        for (int a = 0; a < nAgents; ++a) {
                            int nextAct = Random::nextIntFromDiscrete(actPMF);
                            if (isFermionic && chosenActs[nextAct]) {
                                isValid = false;
                            } else {
                                chosenActs[nextAct] = true;
                                sample[Event<AGENT>(t, agent, nextAct)] = 1;
                                t1State += agent.consequences(nextAct);
                            }
                        }
                        if(!isValid) {agentId = AGENT::domainSize(); t=nTimesteps;}
                    } else {
                        for (int a = 0; a < nAgents; ++a) {
                            int nextAct = Random::nextIntFromDiscrete(actPMF);
                            sample[Event<AGENT>(t, agent, nextAct)] += 1;
                            t1State += agent.consequences(nextAct);
                        }
                    }

//                    ActFermionicDistribution actDistribution(actPMF);
//                    std::vector<bool> chosenActs = actDistribution.sampleUnordered(nAgents);
////                    std::cout << "choosing " << nAgents << " " << chosenActs << std::endl;
//                    if(chosenActs.size() == 0) { // more agents than non-zero probs
//                        isValid = false;
//                        agentId=AGENT::domainSize();
//                        t=nTimesteps;
//                    }

                }
                t0State.setToZero();
                t0State.swap(t1State);
            }
            if (++nAttempts > 40000)
                throw (std::runtime_error(
                        "Can't create act-Fermionic trajectoryPrior sample of Trajectory. Too many agents for Fermionicity to be a good assumption."));
        } while (!isValid);
        return sample;
    }

protected:

    // Pushes the constraint x -> y onto 'constraints'
    // under the assumption that
    // 0 <= x <= 1
    // and 0 <= y_i <= 1
    // by using the identity (see paper)
    //
    void push_xImpliesY(int x, const Constraint<ABM::occupation_type> &y) {
        if(y.upperBound != std::numeric_limits<ABM::occupation_type>::max()) {
            Constraint<ABM::occupation_type> upperBoundConstraint(std::numeric_limits<ABM::occupation_type>::lowest(), 0.0);
            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
                if (y.coefficients.values[i] > 0.0) upperBoundConstraint.upperBound += y.coefficients.values[i];
                upperBoundConstraint.coefficients.insert(y.coefficients.indices[i], y.coefficients.values[i]);
            }
            upperBoundConstraint.coefficients.insert(x, upperBoundConstraint.upperBound - y.upperBound);
            this->addFactor(upperBoundConstraint);
        }
        if(y.lowerBound != std::numeric_limits<ABM::occupation_type>::lowest()) {
            Constraint<ABM::occupation_type> lowerBoundConstraint(std::numeric_limits<ABM::occupation_type>::lowest(), 0.0);
            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
                if (y.coefficients.values[i] < 0.0) lowerBoundConstraint.upperBound -= y.coefficients.values[i];
                lowerBoundConstraint.coefficients.insert(y.coefficients.indices[i], -y.coefficients.values[i]);
            }
            lowerBoundConstraint.coefficients.insert(x, lowerBoundConstraint.upperBound + y.lowerBound);
            this->addFactor(lowerBoundConstraint);
        }
    }

};


#endif //ABMCMC_FORECAST_H
