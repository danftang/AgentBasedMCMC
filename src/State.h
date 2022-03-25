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
public:
    static const std::vector<std::vector<Event<AGENT>>> incomingEventsByState;

    int       time;
    AGENT     agent;

    State(): time(-1), agent(-1) { }

    State(int time, const AGENT &agent):
    time(time),
    agent(agent) { }

    ABM::occupation_type fermionicOccupationUpperBound() const {
        assert(agent < incomingEventsByState.size());
        return AGENT::actDomainSize() < incomingEventsByState[agent].size()?AGENT::actDomainSize():incomingEventsByState[agent].size();
    }


    ABM::occupation_type forwardOccupation(const std::vector<ABM::occupation_type> &trajectory) const {
        ABM::occupation_type occupation = 0;
        int beginIndex = Event<AGENT>(time, agent, 0).id;
        int endIndex = beginIndex + AGENT::actDomainSize();
        for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
            occupation += trajectory[eventId];
        }
        return occupation;
    }


    ABM::occupation_type backwardOccupation(const std::vector<ABM::occupation_type> &trajectory) const {
        ABM::occupation_type occupation = 0;
        for(const Event<AGENT> &incomingEvent : State<AGENT>::incomingEventsByState[agent]) {
            occupation += trajectory[Event<AGENT>(time-1, incomingEvent.agent(), incomingEvent.act()).id];
        }
        return occupation;
    }


    LinearSum<ABM::occupation_type> operator *(ABM::occupation_type c) const {
        LinearSum<ABM::occupation_type> eventVector;
        int beginIndex = Event<AGENT>(time, agent, 0).id;
        for(int actId=0; actId < AGENT::actDomainSize(); ++actId) {
            eventVector += c*X(beginIndex+actId);
        }
        return eventVector;
    }

    friend LinearSum<ABM::occupation_type> operator *(ABM::occupation_type c, const State<AGENT> &state) {
        return state * c;
    }

    friend std::ostream &operator <<(std::ostream &out, const State &state) {
        out << state.agent << "@t" << state.time;
        return out;
    }


//    friend LinearSum &operator +=(LinearSum &sum, const std::pair<double,State<AGENT>> &stateTerm) {
//        int beginIndex = Event(stateTerm.second.time, stateTerm.second.agent, 0);
//        for(int actId=0; actId < AGENT::Act::domainSize; ++actId) {
//            sum += LinearSum(stateTerm.first, beginIndex+actId);
//        }
//        return sum;
//    }


    static std::vector<std::vector<Event<AGENT>>> calculateIncomingEventsByState() {
        std::vector<std::vector<Event<AGENT>>> endStateToEvents(AGENT::domainSize());
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
const std::vector<std::vector<Event<AGENT>>> State<AGENT>::incomingEventsByState = State<AGENT>::calculateIncomingEventsByState();

//template<typename AGENT>
//double State<AGENT>::backwardOccupationNumber(const std::vector<double> &trajectory) const {
//    assert(time != 0);
//    double occupation = 0;
//    for(const Event<AGENT> &incomingEvent : incomingEventsByState[agent]) {
//        occupation += trajectory[Event(time-1, incomingEvent.agent(), incomingEvent.act())];
//    }
//    return occupation;
//}
//
//template<typename AGENT>
//double State<AGENT>::forwardOccupationNumber(const std::vector<double> &trajectory) const {
//    double occupation = 0;
//    for(int act=0; act<AGENT::actDomainSize(); ++act) {
//        occupation += trajectory[Event(time,agent,act)];
//    }
//    return occupation;
//}
//
//template<typename AGENT>
//double State<AGENT>::occupationNumber(const std::vector<double> &trajectory) const {
//    if(trajectory.size() == Trajectory<AGENT>::dimension(time)) return backwardOccupationNumber(trajectory);
//    return forwardOccupationNumber(trajectory);
//}

#endif //GLPKTEST_STATE_H
