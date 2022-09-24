//
// Created by daniel on 23/09/22.
//

#ifndef ABMCMC_EXTENDEDTRAJECTORY2_H
#define ABMCMC_EXTENDEDTRAJECTORY2_H

#include "Trajectory.h"

template<class AGENT, int NTIMESTEPS>
class ExtendedTrajectory2: public Trajectory<AGENT,NTIMESTEPS> {
public:
    typedef typename Trajectory<AGENT,NTIMESTEPS>::occupation_type value_type;
    typedef AGENT agent_type;

    static constexpr size_t dimension = Trajectory<AGENT,NTIMESTEPS>::dimension + NTIMESTEPS*AGENT::domainSize;

    using Trajectory<AGENT,NTIMESTEPS>::operator [];
    using Trajectory<AGENT,NTIMESTEPS>::indexOf;

    ExtendedTrajectory2(): Trajectory<AGENT,NTIMESTEPS>(dimension) { }

protected:
    ExtendedTrajectory2(size_t nElements): Trajectory<AGENT, NTIMESTEPS>(nElements) {}
public:

    const value_type &operator [](const State<AGENT> &state) const {
        return (*this)[indexOf(state)];
    }


//    static SparseVec<value_type> coefficients(const State<AGENT> &state) {
//        SparseVec<value_type> coeffs;
//        coeffs.insert(indexOf(state), 1);
//        return coeffs;
//    }

    static int indexOf(const State<AGENT> &state) {
        return AGENT::domainSize*AGENT::actDomainSize*NTIMESTEPS + state.time*AGENT::domainSize + state.agent;
    }


    static EqualityConstraints<value_type> constraints() {
        EqualityConstraints<value_type> constraints;
        for(int time = 0; time < NTIMESTEPS; ++time) {
            for(int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
                // forward occupation
                SparseVec<ABM::coefficient_type> forwardCoeffs;
                for (int act = 0; act < AGENT::actDomainSize; ++act) {
                    forwardCoeffs.insert(indexOf(Event<AGENT>(time, agentState, act)), 1);
                }
                forwardCoeffs.insert(indexOf(State<agent_type>(time, agentState)), -1);
                constraints.emplace_back(forwardCoeffs, 0);

                // reverse occupation
                if(time > 0) {
                    SparseVec<ABM::coefficient_type> backwardCoeffs;
                    for (const Event<AGENT> &inEdge: State<AGENT>::incomingEventsByState[agentState]) {
                        backwardCoeffs.insert(indexOf(Event<AGENT>(time - 1, inEdge.agent(), inEdge.act())), 1);
                    }
                    backwardCoeffs.insert(indexOf(State<agent_type>(time, agentState)), -1);
                    constraints.emplace_back(backwardCoeffs, 0);
                }
            }
        }
        return constraints;
    }

//    void recalculateDependentVariables() {
//        for(int t=0; t<NTIMESTEPS; ++t) recalculateDependentVariables(t);
//    }

    void recalculateDependentVariables(int timestep) {
        assert(timestep > 0 && timestep < NTIMESTEPS);
        // backward occupation allows us to generate a trajectory one timestep at a time
        for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
            State<AGENT> state(timestep, agentId);
            (*this)[indexOf(state)] = state.backwardOccupation(*this); // stateOccupation;
        }

    }


    ModelState<AGENT> endState() const {
        ModelState<AGENT> endState;
        for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
            endState[agentId] = State<AGENT>(NTIMESTEPS, agentId).backwardOccupation(*this);
        }
        return endState;
    }

    void setStartState(const ModelState<AGENT> &startState) {
        for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
            (*this)[indexOf(State<AGENT>(0,agentId))] = startState[agentId];
        }
    }
};


#endif //ABMCMC_EXTENDEDTRAJECTORY2_H
