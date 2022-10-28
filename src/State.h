//
// Created by daniel on 04/05/2021.
//

#ifndef GLPKTEST_STATE_H
#define GLPKTEST_STATE_H

#include <utility>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <assert.h>
#include "LinearSum.h"
#include "Event.h"
#include "ABM.h"

template<typename AGENT>
class State {
protected:
    static const std::vector<std::vector<Event<AGENT>>> incomingEventsAtT0ByAgentId; // evants at time t=0
public:

    int       time;
    AGENT     agent;

    State(): time(-1), agent(-1) { }

    State(int time, const AGENT &agent):
    time(time),
    agent(agent) { }

    explicit State(int stateTrajectoryIndex): State(std::div(stateTrajectoryIndex, AGENT::domainSize)) {}

    explicit State(std::div_t agentIdDivDomainSize): time(agentIdDivDomainSize.quot), agent(agentIdDivDomainSize.rem) { }


    int id() const {
        return time*AGENT::domainSize + agent;
    }


//    ABM::occupation_type fermionicOccupationUpperBound() const {
//        assert(agent < incomingEventsByState.size());
//        return AGENT::actDomainSize < incomingEventsByState[agent].size()?AGENT::actDomainSize:incomingEventsByState[agent].size();
//    }

    template<class DOMAIN>
    typename DOMAIN::value_type forwardOccupation(const DOMAIN &trajectory) const {
        typename DOMAIN::value_type occupation = 0;
        for(int actId = 0; actId < AGENT::actDomainSize; ++actId) {
            occupation += trajectory[DOMAIN::indexOf(Event<AGENT>(time, agent, actId))];
        }
//        int beginIndex = Event<AGENT>(time, agent, 0).id;
//        int endIndex = beginIndex + AGENT::actDomainSize;
//        for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
//            occupation += trajectory[eventId];
//        }
        return occupation;
    }


    template<class DOMAIN>
    typename DOMAIN::value_type backwardOccupation(const DOMAIN &trajectory) const {
        typename DOMAIN::value_type occupation = 0;
        for(const Event<AGENT> &incomingEvent : incomingEventsAtT0ByAgentId[agent]) {
            occupation += trajectory[DOMAIN::indexOf(Event<AGENT>(time-1, incomingEvent.agent(), incomingEvent.act()))];
        }
        return occupation;
    }

    std::vector<Event<AGENT>> forwardOccupationDependencies() const {
        std::vector<Event<AGENT>> dependencies;
//        int beginIndex = Event<AGENT>(time, agent, 0).id;
//        int endIndex = beginIndex + AGENT::actDomainSize;
//        for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
//            dependencies.push_back(eventId);
//        }
        for(int actId=0; actId < AGENT::actDomainSize; ++actId) {
            dependencies.push_back(Event<AGENT>(time,agent,typename AGENT::Act(actId)));
        }
        return dependencies;
    }

    std::vector<Event<AGENT>> backwardOccupationDependencies() const {
        assert(time > 0);
        std::vector<Event<AGENT>> dependencies;
        for(const Event<AGENT> &t0Event: incomingEventsAtT0ByAgentId[agent]) {
            dependencies.push_back(Event<AGENT>(time-1, t0Event.agent(), t0Event.act()));
        }
        return dependencies;
    }


//    ABM::occupation_type fermionicBoundedForwardOccupation(const std::vector<ABM::occupation_type> &trajectory) const {
//        ABM::occupation_type occupation = 0;
//        int beginIndex = Event<AGENT>(time, agent, 0).id;
//        int endIndex = beginIndex + AGENT::actDomainSize;
//        for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
//            occupation += trajectory[eventId]>0?1:0;
//        }
//        return occupation;
//    }


//    LinearSum<ABM::occupation_type> operator *(ABM::occupation_type c) const {
//        LinearSum<ABM::occupation_type> eventVector;
//        int beginIndex = Event<AGENT>(time, agent, 0).id;
//        for(int actId=0; actId < AGENT::actDomainSize; ++actId) {
//            eventVector += c*X(beginIndex+actId);
//        }
//        return eventVector;
//    }

//    friend LinearSum<ABM::occupation_type> operator *(ABM::occupation_type c, const State<AGENT> &state) {
//        return state * c;
//    }

    friend std::ostream &operator <<(std::ostream &out, const State &state) {
        out << state.agent << "@t" << state.time;
        return out;
    }


    static std::vector<std::vector<Event<AGENT>>> calculateIncomingEventsByAgentId() {
        std::vector<std::vector<Event<AGENT>>> consequenceToEvent(AGENT::domainSize);
        for(int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
            AGENT agent(agentState);
            for (int act = 0; act < AGENT::actDomainSize; ++act) {
                for(const AGENT &endState: agent.consequences(act)) {
                    consequenceToEvent[endState].push_back(Event(0, agent, act));
                }
            }
        }
        return consequenceToEvent;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void load(Archive &ar, const unsigned int version) {
        int agentId;
        ar & time & agentId;
        agent = AGENT(agentId);
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        int agentId = agent;
        ar & time & agentId;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();


};

template<typename AGENT>
const std::vector<std::vector<Event<AGENT>>> State<AGENT>::incomingEventsAtT0ByAgentId = State<AGENT>::calculateIncomingEventsByAgentId();

#endif //GLPKTEST_STATE_H
