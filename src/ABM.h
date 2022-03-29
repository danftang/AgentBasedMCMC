// Static class that contains functions that create ConvexPolyhedrons for the
// various constarints that must be satisfied by valid ABM trajectories.
//
// Created by daniel on 02/08/2021.
//

#ifndef GLPKTEST_ABMCONSTRAINTS_H
#define GLPKTEST_ABMCONSTRAINTS_H

//#include "ConvexPolyhedron.h"
//#include "Event.h"

class ABM {
public:
    typedef int occupation_type;





//    // All constraints for an act fermionic trajectory
//    template<typename AGENT>
//    static ConvexPolyhedron<occupation_type> actFermionicABMConstraints(int nTimesteps) {
//        return  continuityConstraints<AGENT>(nTimesteps) +
//                interactionConstraints<AGENT>(nTimesteps) +
//                actFermionicConstraints<AGENT>(nTimesteps);
//    }
//
//
//    template<typename AGENT>
//    static ConvexPolyhedron<occupation_type> continuityConstraints(int nTimesteps) {
//        ConvexPolyhedron<occupation_type> constraints;
//        std::vector<std::vector<Event<AGENT>>> incomingEdges = consequencesByEndState<AGENT>();
//        for(int time = 1; time < nTimesteps; ++time) {
//            for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
//                Constraint<occupation_type> &constraint = constraints.emplace_back(0.0 ,0.0);
//                // outgoing edges
//                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
//                    constraint.coefficients.insert(Event<AGENT>(time, agentState, act), 1.0);
//                }
//                // incoming edges
////                int timeOffset = (time-1)*AGENT::domainSize()*AGENT::actDomainSize();
//                for (const Event<AGENT> &inEdge: incomingEdges[agentState]) {
//                    constraint.coefficients.insert(Event<AGENT>(time-1,inEdge.agent(),inEdge.act()), -1.0);
//                }
//            }
//        }
//        return constraints;
//    }
//
//
//    template<typename AGENT>
//    static ConvexPolyhedron<occupation_type> interactionConstraints(int nTimesteps) {
//        ConvexPolyhedron<occupation_type> constraints;
//        for(int time = 0; time < nTimesteps; ++time) {
//            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
//                AGENT agent(agentState);
//                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
//                    for(const Constraint<occupation_type> &actConstraint : agent.constraints(time, act)) {
//                        push_xImpliesY(constraints, Event(time,agent,act), actConstraint);
//                    }
//                }
//            }
//        }
//        return constraints;
//    }
//
//
//    template<typename AGENT>
//    static ConvexPolyhedron<occupation_type> actFermionicConstraints(int nTimesteps) {
//        ConvexPolyhedron<occupation_type> constraints;
//        for(int time = 0; time < nTimesteps; ++time) {
//            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
//                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
//                    Event<AGENT> event = Event<AGENT>(time,agentState,act);
//                    constraints.push_back(0.0 <= 1.0*event <= 1.0);
//                }
//            }
//        }
//        return constraints;
//    }
//
//
//    template<typename AGENT>
//    static ConvexPolyhedron<occupation_type> stateOccupationNumbersAsAuxiliary(int nTimesteps) {
//        ConvexPolyhedron<occupation_type> constraints;
//        for(int time = 0; time < nTimesteps; ++time) {
//            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
//                Constraint<occupation_type> & occupation = constraints.emplace_back(0.0, AGENT::actDomainSize()); // assumes act Fermionicity
//                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
//                    occupation.coefficients.insert(Event<AGENT>(time,agentState,act), 1.0);
//                }
//            }
//        }
//        return constraints;
//    }
//
//
////    static ConvexPolyhedron<occupation_type> startStateConstraintsToTrajectoryConstraints(const ConvexPolyhedron<occupation_type> &startStateConstraints) {
////        ConvexPolyhedron<occupation_type> trajectoryConstraints;
////        for(const Constraint<occupation_type> &constraint: startStateConstraints) {
////            Constraint<occupation_type> trajectoryConstraint = startStateConstraintToTrajectoryConstraint(constraint);
////            if(trajectoryConstraint.coefficients.sparseSize() != 0) {
////                trajectoryConstraints.push_back(std::move(trajectoryConstraint));
////            }
////        }
////        return trajectoryConstraints;
////    }
//
//protected:
//
//    // Pushes the constraint x -> y onto 'constraints'
//    // under the assumption that
//    // 0 <= x <= 1
//    // and 0 <= y_i <= 1
//    // by using the identity (see paper)
//    //
//    static void push_xImpliesY(std::vector<Constraint<occupation_type>> &constraints, int x, const Constraint<occupation_type> &y) {
//        if(y.upperBound != INFINITY) {
//            Constraint<occupation_type> &upperBoundConstraint = constraints.emplace_back(-INFINITY, 0.0);
//            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
//                if (y.coefficients.values[i] > 0.0) upperBoundConstraint.upperBound += y.coefficients.values[i];
//                upperBoundConstraint.coefficients.insert(y.coefficients.indices[i], y.coefficients.values[i]);
//            }
//            upperBoundConstraint.coefficients.insert(x, upperBoundConstraint.upperBound - y.upperBound);
//        }
//        if(y.lowerBound != -INFINITY) {
//            Constraint<occupation_type> &lowerBoundConstraint = constraints.emplace_back(-INFINITY, 0.0);
//            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
//                if (y.coefficients.values[i] < 0.0) lowerBoundConstraint.upperBound -= y.coefficients.values[i];
//                lowerBoundConstraint.coefficients.insert(y.coefficients.indices[i], -y.coefficients.values[i]);
//            }
//            lowerBoundConstraint.coefficients.insert(x, lowerBoundConstraint.upperBound + y.lowerBound);
//        }
//    }
//
//
//    template<typename AGENT>
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
//
//
////    template<typename AGENT>
////    static Constraint<occupation_type> startStateConstraintToTrajectoryConstraint(const Constraint<occupation_type> &startStateConstraint) {
////        LinearSum<occupation_type> trajectoryCoeffs;
////        double fermionicUpperBound = 0;
////        for(int i=0; i<startStateConstraint.coefficients.sparseSize(); ++i) {
////            State<AGENT> state(0, startStateConstraint.coefficients.indices[i]);
////            trajectoryCoeffs += startStateConstraint.coefficients.values[i]*state;
////            fermionicUpperBound += state.occupationUpperBound();
////        }
////        if(startStateConstraint.lowerBound > 0.0 || startStateConstraint.upperBound < fermionicUpperBound) {
////            return std::max(startStateConstraint.lowerBound,0.0) <= trajectoryCoeffs <= std::min(startStateConstraint.upperBound,fermionicUpperBound);
////        }
////        return Constraint();
////    }

};


#endif //GLPKTEST_ABMCONSTRAINTS_H
