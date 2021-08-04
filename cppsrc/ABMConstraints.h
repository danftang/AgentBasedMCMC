//
// Created by daniel on 02/08/2021.
//

#ifndef GLPKTEST_ABMCONSTRAINTS_H
#define GLPKTEST_ABMCONSTRAINTS_H

template<typename AGENT>
class ABMConstraints {
public:

    static ConvexPolyhedron actFermionicABMConstraints(int nTimesteps) {
        return  continuityConstraints(nTimesteps) +
                interactionConstraints(nTimesteps) +
                actFermionicConstraints(nTimesteps);
    }

    static ConvexPolyhedron continuityConstraints(int nTimesteps) {
        ConvexPolyhedron constraints;
        std::vector<std::vector<int>> incomingEdges = consequencesByEndState();
        for(int time = 1; time < nTimesteps; ++time) {
            for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                glp::Constraint &constraint = constraints.emplace_back(0.0 ,0.0);
                // outgoing edges
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    constraint.coefficients.insert(Event<AGENT>(time, agentState, act), 1.0);
                }
                // incoming edges
                int timeOffset = (time-1)*AGENT::domainSize()*AGENT::actDomainSize();
                for (int inEdge: incomingEdges[agentState]) {
                    constraint.coefficients.insert(timeOffset + inEdge, -1.0);
                }
            }
        }
        return constraints;
    }

    static ConvexPolyhedron interactionConstraints(int nTimesteps) {
        ConvexPolyhedron constraints;
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                AGENT agent(agentState);
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    for(const glp::Constraint &actConstraint : agent.constraints(time, act)) {
                        push_xImpliesY(constraints, Event(time,agent,act), actConstraint);
                    }
                }
            }
        }
        return constraints;
    }

    // Returns constraint x -> y
    // under the assumption that
    // 0 <= x <= 1
    // and 0 <= y_i <= 1
    // by using the identity
    //
    static void push_xImpliesY(std::vector<glp::Constraint> &constraints, int x, const glp::Constraint &y) {
        if(y.upperBound != INFINITY) {
            glp::Constraint &upperBoundConstraint = constraints.emplace_back(-INFINITY, 0.0);
            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
                if (y.coefficients.values[i] > 0.0) upperBoundConstraint.upperBound += y.coefficients.values[i];
                upperBoundConstraint.coefficients.insert(y.coefficients.indices[i], y.coefficients.values[i]);
            }
            upperBoundConstraint.coefficients.insert(x, upperBoundConstraint.upperBound - y.upperBound);
        }
        if(y.lowerBound != -INFINITY) {
            glp::Constraint &lowerBoundConstraint = constraints.emplace_back(-INFINITY, 0.0);
            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
                if (y.coefficients.values[i] < 0.0) lowerBoundConstraint.upperBound -= y.coefficients.values[i];
                lowerBoundConstraint.coefficients.insert(y.coefficients.indices[i], -y.coefficients.values[i]);
            }
            lowerBoundConstraint.coefficients.insert(x, lowerBoundConstraint.upperBound + y.lowerBound);
        }
    }


    static ConvexPolyhedron actFermionicConstraints(int nTimesteps) {
        ConvexPolyhedron constraints;
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    Event<AGENT> event = Event<AGENT>(time,agentState,act);
                    constraints.push_back(0.0 <= 1.0*event <= 1.0);
                }
            }
        }
        return constraints;
    }


    static std::vector<std::vector<int>> consequencesByEndState() {
        std::vector<std::vector<int> > endStateToEvents(AGENT::domainSize());
        std::vector<AGENT> consequences;
        for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
            AGENT agent(agentState);
            for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                consequences = agent.consequences(act);
                for(const AGENT &endState: consequences) {
                    endStateToEvents[endState].push_back(Event(0,agent,act));
                }
            }
        }
        return endStateToEvents;
    }


};


#endif //GLPKTEST_ABMCONSTRAINTS_H
