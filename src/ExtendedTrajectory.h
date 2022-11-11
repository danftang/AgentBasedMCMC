//
// Created by daniel on 23/09/22.
//

#ifndef ABMCMC_EXTENDEDTRAJECTORY_H
#define ABMCMC_EXTENDEDTRAJECTORY_H

#include "Trajectory.h"

template<class AGENT, int NTIMESTEPS>
class ExtendedTrajectory: public Trajectory<AGENT,NTIMESTEPS> {
public:
    typedef typename Trajectory<AGENT,NTIMESTEPS>::value_type value_type;
    typedef AGENT agent_type;

    static constexpr size_t dimension = Trajectory<AGENT,NTIMESTEPS>::dimension + NTIMESTEPS*AGENT::domainSize;

    using Trajectory<AGENT,NTIMESTEPS>::operator [];
    using Trajectory<AGENT,NTIMESTEPS>::indexOf;

    ExtendedTrajectory(): Trajectory<AGENT,NTIMESTEPS>(dimension) { }
    ExtendedTrajectory(const SparseVec<int> &sparseInit): ExtendedTrajectory(dimension, sparseInit) { }


protected:
    ExtendedTrajectory(size_t nElements): Trajectory<AGENT, NTIMESTEPS>(nElements) {}
    ExtendedTrajectory(size_t nElements, const SparseVec<int> &sparseInit): Trajectory<AGENT,NTIMESTEPS>(dimension, sparseInit) { }
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
            for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
                State<AGENT> state(time,agentId);

                if(time == 0) {
                    SparseVec<value_type> forwardCoeffs;
                    for (int act = 0; act < AGENT::actDomainSize; ++act) {
                        forwardCoeffs.insert(indexOf(Event<AGENT>(time, agentId, act)), 1);
                    }
                    forwardCoeffs.insert(indexOf(state), -1);
                    constraints.emplace_back(forwardCoeffs, 0);
                } else {
                    SparseVec<value_type> backwardCoeffs;
                    for (const Event<AGENT> &inEdge: state.backwardOccupationDependencies()) {
                        backwardCoeffs.insert(indexOf(inEdge),1);
                    }
                    backwardCoeffs.insert(indexOf(state), -1);
                    constraints.emplace_back(backwardCoeffs, 0);

                    SparseVec<value_type> doubleCoeffs;
                    for (const Event<AGENT> &inEdge: state.backwardOccupationDependencies()) {
                        doubleCoeffs.insert(indexOf(inEdge),1);
                    }
                    for (int act = 0; act < AGENT::actDomainSize; ++act) {
                        doubleCoeffs.insert(indexOf(Event<AGENT>(time, agentId, act)), -1);
                    }
                    constraints.emplace_back(doubleCoeffs, 0);
                }

//                // forward occupation
//                SparseVec<value_type> forwardCoeffs;
//                for (int act = 0; act < AGENT::actDomainSize; ++act) {
//                   forwardCoeffs.insert(indexOf(Event<AGENT>(time, agentId, act)), 1);
//                }
//                forwardCoeffs.insert(indexOf(state), -1);
//                constraints.emplace_back(forwardCoeffs, 0);
//                // reverse occupation
//                if(time > 0) {
//                    SparseVec<value_type> backwardCoeffs;
//                    for (const Event<AGENT> &inEdge: state.backwardOccupationDependencies()) {
//                        backwardCoeffs.insert(indexOf(inEdge),1);
//                    }
//                    backwardCoeffs.insert(indexOf(state), -1);
//                    constraints.emplace_back(backwardCoeffs, 0);
//                }

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

    void sanityCheck() const {
        for(int t=0; t<NTIMESTEPS; ++t) {
            for(int agentId=0; agentId < AGENT::domainSize; ++agentId) {
                int stateOccupation = (*this)[State<AGENT>(t,agentId)];
                int sumOfEvents = 0;
                for(int actId = 0; actId < AGENT::actDomainSize; ++actId) {
                    sumOfEvents += (*this)[Event<AGENT>(t,agentId,actId)];
                }
                assert(sumOfEvents == stateOccupation);
            }
        }
    }
};


#endif //ABMCMC_EXTENDEDTRAJECTORY_H
