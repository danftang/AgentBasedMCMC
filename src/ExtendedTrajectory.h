//
// Created by daniel on 18/09/22.
//

#ifndef ABMCMC_EXTENDEDTRAJECTORY_H
#define ABMCMC_EXTENDEDTRAJECTORY_H

#include "Trajectory.h"

template<class AGENT>
class ExtendedTrajectory {
public:
    typedef ABM::occupation_type    value_type;
    typedef AGENT                   agent_type;

    Trajectory<AGENT> actTrajectory;
    std::vector<ModelState<AGENT>> stateTrajectory; // not the end state

    class TrajectoryAccessor {
    public:
        typedef typename Trajectory<AGENT>::value_type value_type;

        value_type &actEntry;
        value_type &modelStateEntry;

        TrajectoryAccessor(value_type &actEntry, value_type &modelStateEntry):
        actEntry(actEntry), modelStateEntry(modelStateEntry) {}

        operator const value_type &() const {
            return actEntry;
        }

        TrajectoryAccessor &operator =(const value_type &newValue) {
            modelStateEntry += newValue - actEntry;
            actEntry = newValue;
            return *this;
        }

        TrajectoryAccessor &operator +=(const value_type &increment) {
            modelStateEntry += increment;
            actEntry += increment;
            return *this;
        }

    };


    ExtendedTrajectory(int nTimesteps): actTrajectory(nTimesteps), stateTrajectory(nTimesteps) { }


    size_t size() const { return nTimesteps()*(AGENT::domainSize*(AGENT::actDomainSize+1));  }

    size_t nTimesteps() const { return stateTrajectory.size(); }

    ModelState<AGENT> modelState(int time) const { return ModelState<AGENT>(*this, time); }

    // model state access (const only as modification of state is ambiguous)
    value_type operator [](const State<AGENT> &state) const {
        if(state.time == nTimesteps()) return actTrajectory[state]; // we don't store the end state
        return stateTrajectory[state.time][state.agent];
    }

    // trajectory access
    const value_type &operator [](const Event<AGENT> &event) const {
        return actTrajectory[event];
    };

    TrajectoryAccessor operator [](const Event<AGENT> &event) {
        return TrajectoryAccessor(actTrajectory[event], stateTrajectory[event.time()][event.agent()]);
    };

    // direct element access
    const value_type &operator [](int index) const {
        div_t div = std::div(index, AGENT::domainSize*(AGENT::actDomainSize+1));
        int time = div.quot;
        if(div.rem < AGENT::domainSize) {
            return stateTrajectory[time][div.rem];
        }
        return actTrajectory[AGENT::domainSize * AGENT::actDomainSize * time + div.rem - AGENT::domainSize];
    }


    value_type &operator [](int index) {
        div_t div = std::div(index, AGENT::domainSize*(AGENT::actDomainSize+1));
        int time = div.quot;
        if(div.rem < AGENT::domainSize) {
            return stateTrajectory[time][div.rem];
        }
        return actTrajectory[AGENT::domainSize * AGENT::actDomainSize * time + div.rem - AGENT::domainSize];
    }

    static SparseVec<value_type> coefficients(const State<AGENT> &state) {
        SparseVec<value_type> coeffs;
        coeffs.insert(indexOf(state), 1);
        return coeffs;
    }

    static int indexOf(const Event<AGENT> &event) {
        return event.time() * (AGENT::domainSize*(AGENT::actDomainSize+1)) + AGENT::domainSize + event.agent()*AGENT::actDomainSize + event.act();
    }

    static int indexOf(const State<AGENT> &state) {
        return state.time * (AGENT::domainSize*(AGENT::actDomainSize+1)) + state.agent;
    }


    static EqualityConstraints<value_type> constraints(int nTimesteps) {
        EqualityConstraints<value_type> constraints;
        for(int time = 0; time < nTimesteps; ++time) {
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


};


#endif //ABMCMC_EXTENDEDTRAJECTORY_H
