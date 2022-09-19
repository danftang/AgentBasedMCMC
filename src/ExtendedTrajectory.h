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

    Trajectory<AGENT> trajectory;
    std::vector<ModelState<AGENT>> modelState;

    class TrajectoryAccessor {
    public:
        typedef typename Trajectory<AGENT>::value_type value_type;

        value_type &actEntry;
        value_type &modelStateEntry;

        TrajectoryAccessor(value_type &actEntry, value_type &modelStateEntry):
        actEntry(actEntry), modelStateEntry(modelStateEntry) {}

        operator const value_type &() {
            return actEntry;
        }

        value_type &operator =(const value_type &newValue) {
            modelStateEntry += newValue - actEntry;
            actEntry = newValue;
        }
    };


    // model state access (const only as modification of state is ambiguous)
    const value_type &operator [](const State<AGENT> &state) const {
        return modelState[state.time][state.agent];
    }

    // trajectory access
    const value_type &operator [](const Event<AGENT> &event) const {
        return trajectory[event];
    };

    TrajectoryAccessor operator [](const Event<AGENT> &event) {
        return TrajectoryAccessor(trajectory[event], modelState[event.time][event.agent()]);
    };

    // direct element access
    const value_type &operator [](int index) const {
        if(index < trajectory.size()) return trajectory[index];
        int stateTrajectoryIndex = index - trajectory.size();
        int time = stateTrajectoryIndex / AGENT::domainSize();
        int state = stateTrajectoryIndex % AGENT::domainSize();
        return modelState[time][state];
    }


    value_type &operator [](int index) {
        if(index < trajectory.size()) return trajectory[index];
        int stateTrajectoryIndex = index - trajectory.size();
        int time = stateTrajectoryIndex / AGENT::domainSize();
        int state = stateTrajectoryIndex % AGENT::domainSize();
        return modelState[time][state];
    }


};


#endif //ABMCMC_EXTENDEDTRAJECTORY_H
