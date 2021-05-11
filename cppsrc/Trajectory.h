//
// Created by daniel on 04/05/2021.
//

#ifndef GLPKTEST_TRAJECTORY_H
#define GLPKTEST_TRAJECTORY_H

#include "glpkpp.h"

template<typename AGENT>
class Trajectory: public glp::SparseVec {
public:
    Trajectory() { }
    Trajectory(glp::SparseVec &&rvalue): glp::SparseVec(0) {
        swap(rvalue);
    }

    Trajectory(const glp::SparseVec &lvalue): glp::SparseVec(lvalue) { }

    Trajectory &operator =(glp::SparseVec &&rvalue) {
        swap(rvalue);
        return *this;
    }

    Trajectory &operator =(const glp::SparseVec &lvalue) {
        glp::SparseVec::operator=(lvalue);
        return *this;
    }

    // event count
    double operator()(int time, const AGENT &agent, const typename AGENT::Act &act) const {
        return (*this)[Event(time,agent,act)];
    }

    // occupation number
    double operator()(int time, const AGENT &agent) const {
        int beginIndex = Event(time,agent,0);
        int endIndex = beginIndex + AGENT::actDomainSize();
        double occupation = 0.0;
        int eventId;
        for(int i=1; i<=sparseSize(); ++i) {
            eventId = indices[i];
            if(beginIndex <= eventId && eventId < endIndex) occupation += values[i];
        }
        return occupation;
    }

    std::vector<std::map<AGENT,double>> getStateTrajectory(int time) {
        std::vector<std::map<AGENT,double>> stateTrajectory;
        int eventId;
        for(int i=1; i<=sparseSize(); ++i) {
            auto event = Event<AGENT>(indices[i]);
            if(event.time() >= stateTrajectory.size()) stateTrajectory.resize(event.time()+1);
            stateTrajectory[event.time()][event.agent()] += values[i];
        }
        return stateTrajectory;
    }

    double logPrior() {
        double logP = 0.0;
        std::vector<std::map<AGENT,double>> stateTrajectory = getStateTrajectory();
        for(auto eventOccupation : *this) {
            auto event = Event<AGENT>(eventOccupation.index);
            logP += eventOccupation.value * ln( event.agent().timestep(stateTrajectory[event.time()])[event.act()] );
//                   - CombinatoricsUtils.factorialLog(occupation) // TODO: add this for non-state-fermionic trajectories
        }
// TODO: Add this as we're no longer assuming state-fermionicity
//        for(state in stateTrajectory) {
//            for((_, occupation) in state.entries) {
//                logP += CombinatoricsUtils.factorialLog(occupation)
//            }
//        }
        return logP;
    }


    friend std::ostream &operator <<(std::ostream &out, const Trajectory<AGENT> &trajectory) {
        Trajectory sortedTrajectory(trajectory);
        sortedTrajectory.sort();

        int timestep = -1;
        for(auto eventCount: sortedTrajectory) {
            Event<AGENT> event = eventCount.index;
            if(event.time() != timestep) {
                out << "time = " << event.time() << std::endl;
                timestep = event.time();
            }
            out << "  " << eventCount.value << "*" << event << std::endl;
        }

        return out;
    }

//    val stateTrajectory: List<Multiset<AGENT>> by lazy {
//            val stateTrajectory = ArrayList<Multiset<AGENT>>()
//            for((event, occupation) in events) {
//                while(stateTrajectory.size <= event.time) stateTrajectory.add(Multiset())
//                stateTrajectory[event.time][event.agent] += occupation
//            }
//            stateTrajectory
//    }
//
//    val finalState: Multiset<AGENT> by lazy {
//            val finalState = Multiset<AGENT>()
//            val penultimateTime = stateTrajectory.lastIndex
//            stateTrajectory
//            .lastOrNull()
//            ?.members
//            ?.forEach { agent ->
//                for(act in model.actDomain) {
//                    val occupation = this[penultimateTime, agent, act]
//                    if(occupation > 0) {
//                        finalState += model.consequences(agent, act) * occupation
//                    }
//                }
//            }
//            finalState
//    }
//
//
//    val events: Sequence<AbstractMap.SimpleEntry<ABM.Event<AGENT, ACT>, Int>>
//    get() = eventVector.nonZeroEntries
//            .asSequence()
//            .map { AbstractMap.SimpleEntry(model.eventDomain[it.key], it.value.toDouble().roundToInt()) }
//
//
//    inline operator fun get(time: Int, agent: AGENT, act: ACT): Int {
//        return get(ABM.Event(time, agent, act))
//    }
//
//    inline operator fun get(event: ABM.Event<AGENT, ACT>): Int {
//        return eventVector[event.ordinal].toDouble().roundToInt()
//    }
//
//
//
//    fun stateAt(time: Int): Multiset<AGENT> {
//        if(time < stateTrajectory.size) return stateTrajectory[time]
//        if(time == stateTrajectory.size && time != 0) return finalState
//        return emptyMultiset()
//    }
//
//
//    fun nAgents(time: Int, agent: AGENT): Int {
//        if(time < stateTrajectory.size) return stateTrajectory[time][agent]
//        if(time == stateTrajectory.size && time != 0) return finalState[agent]
//        return 0
//    }
//
//

};


#endif //GLPKTEST_TRAJECTORY_H
