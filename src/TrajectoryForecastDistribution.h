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

#ifndef ABMCMC_TRAJECTORYFORECASTDISTRIBUTION_H
#define ABMCMC_TRAJECTORYFORECASTDISTRIBUTION_H


#include "WeightedFactoredConvexDistribution.h"
#include "TrajectoryImportance.h"
#include "ABM.h"

template <class AGENT>
class TrajectoryForecastDistribution: public WeightedFactoredConvexDistribution<ABM::occupation_type> {
public:

    explicit TrajectoryForecastDistribution(int nTimesteps):
    WeightedFactoredConvexDistribution<ABM::occupation_type>( [nTimesteps]() {
        return std::unique_ptr<PerturbableFunction<ABM::occupation_type,double>>(new TrajectoryImportance<AGENT>(nTimesteps));
    }) {
        addActFermionicFactors(nTimesteps);
        addContinuityConstraints(nTimesteps);
        addInteractionConstraints(nTimesteps);
    }


    void addContinuityConstraints(int nTimesteps) {
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


    void addInteractionConstraints(int nTimesteps) {
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


    void addActFermionicFactors(int nTimesteps) {
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


//    static ConvexPolyhedron<ABM::occupation_type> stateOccupationNumbersAsAuxiliary(int nTimesteps) {
//        ConvexPolyhedron<ABM::occupation_type> constraints;
//        for(int time = 0; time < nTimesteps; ++time) {
//            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
//                Constraint<ABM::occupation_type> & occupation = constraints.emplace_back(0.0, AGENT::actDomainSize()); // assumes act Fermionicity
//                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
//                    occupation.coefficients.insert(Event<AGENT>(time,agentState,act), 1.0);
//                }
//            }
//        }
//        return constraints;
//    }



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


//    static std::vector<std::vector<Event<AGENT>>> consequencesByEndState() {
//        std::vector<std::vector<Event<AGENT>>> endStateToEvents(AGENT::domainSize());
//        std::vector<AGENT> consequences;
//        for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
//            AGENT agent(agentState);
//            for (int act = 0; act < AGENT::actDomainSize(); ++act) {
//                consequences = agent.consequences(act);
//                for(const AGENT &endState: consequences) {
//                    endStateToEvents[endState].push_back(Event(0,agent,act));
//                }
//            }
//        }
//        return endStateToEvents;
//    }

};


#endif //ABMCMC_TRAJECTORYFORECASTDISTRIBUTION_H
