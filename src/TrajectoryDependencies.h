//
// Created by daniel on 22/09/22.
//

#ifndef ABMCMC_TRAJECTORYDEPENDENCIES_H
#define ABMCMC_TRAJECTORYDEPENDENCIES_H

#include <vector>
#include "Event.h"
#include "State.h"

template<class AGENT>
class TrajectoryDependencies {
public:

    std::vector<State<AGENT>> states;
    std::vector<Event<AGENT>> events;

    template<class STATES = std::initializer_list<State<AGENT>>, class EVENTS = std::initializer_list<Event<AGENT>>>
    TrajectoryDependencies(STATES &&states, EVENTS &&events): states(std::forward<STATES>(states)), events(std::forward<EVENTS>(events)) {
    }

//    template<class TRAJECTORY>
//    std::vector<int> indexDependencies() {
//        std::vector<int> indexDependencies;
//        for(const State<AGENT> &state: states) {
//            for(int index: TRAJECTORY::coefficients(state).indices) indexDependencies.push_back(index);
//        }
//        for(const Event<AGENT> &event: events) indexDependencies.push_back(TRAJECTORY::indexOf(event));
//        return indexDependencies;
//    }
};


#endif //ABMCMC_TRAJECTORYDEPENDENCIES_H
