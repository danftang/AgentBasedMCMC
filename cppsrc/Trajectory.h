//
// Created by daniel on 04/05/2021.
//

#ifndef GLPKTEST_TRAJECTORY_H
#define GLPKTEST_TRAJECTORY_H

#include "glpkpp.h"
#include "ModelState.h"
#include "StateTrajectory.h"
#include "Random.h"

template<typename AGENT>
class Trajectory: public std::vector<double> {
public:
    Trajectory(int nTimesteps): std::vector<double>(nTimesteps*AGENT::domainSize()*AGENT::actDomainSize()+1) { }
    Trajectory(std::vector<double> &&rvalue): std::vector<double>(std::move(rvalue)) { }
    Trajectory(const std::vector<double> &lvalue): std::vector<double>(lvalue) { }

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
        for(int eventId=beginIndex; eventId<endIndex; ++eventId) {
            occupation += (*this)[eventId];
        }
        return occupation;
    }


//    double logPrior() {
//        double logP = 0.0;
//        std::vector<std::map<AGENT,double>> stateTrajectory = toStateTrajectory();
//        for(auto eventOccupation : *this) {
//            auto event = Event<AGENT>(eventOccupation.index);
//            logP += eventOccupation.value * ln( event.agent().timestep(stateTrajectory[event.time()])[event.act()] );
////                   - CombinatoricsUtils.factorialLog(occupation) // TODO: add this for non-state-fermionic trajectories
//        }
//// TODO: Add this as we're no longer assuming state-fermionicity
////        for(state in stateTrajectory) {
////            for((_, occupation) in state.entries) {
////                logP += CombinatoricsUtils.factorialLog(occupation)
////            }
////        }
//        return logP;
//    }

    static Trajectory<AGENT> run(const ModelState<AGENT> &startState, int nTimesteps);


//    friend std::ostream &operator <<(std::ostream &out, const Trajectory<AGENT> &trajectory) {
//        Trajectory sortedTrajectory(trajectory);
//        sortedTrajectory.sort();
//
//        int timestep = -1;
//        for(auto eventCount: sortedTrajectory) {
//            Event<AGENT> event = eventCount.index;
//            if(event.time() != timestep) {
//                out << "time = " << event.time() << std::endl;
//                timestep = event.time();
//            }
//            out << "  " << eventCount.value << "*" << event << std::endl;
//        }
//        return out;
//    }


};

template<typename AGENT>
Trajectory<AGENT> Trajectory<AGENT>::run(const ModelState<AGENT> &startState, int nTimesteps) {
    ModelState<AGENT> t0State = startState;
    ModelState<AGENT> t1State;
    Trajectory<AGENT> trajectory(nTimesteps);
    for(int t=0; t<nTimesteps; ++t) {
        for(const auto &[agent, occupation]: t0State) {
            std::vector<double> actPMF = agent.timestep(t0State, 0.0);
            assert(occupation <= actPMF.size());
            for(int nthAgent=0; nthAgent < occupation; ++nthAgent) {
                typename AGENT::Act act = Random::choose(actPMF);
                trajectory[Event(t,agent,act)] = 1.0;
                t1State += agent.consequences(act);
                actPMF[act] = 0.0;
            }
        }
        t0State.clear();
        t0State.swap(t1State);
    }
    return trajectory;
}


#endif //GLPKTEST_TRAJECTORY_H
