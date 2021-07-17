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
//    static constexpr double tol = SimplexMCMC::tol;

    explicit Trajectory(int nTimesteps): std::vector<double>(nTimesteps*AGENT::domainSize()*AGENT::actDomainSize()+1) { }
    explicit Trajectory(std::vector<double> &&rvalue): std::vector<double>(std::move(rvalue)) { }
    explicit Trajectory(const std::vector<double> &lvalue): std::vector<double>(lvalue) { }

    // event count
    double operator()(int time, const AGENT &agent, const typename AGENT::Act &act) const {
        return (*this)[Event(time,agent,act)];
    }

    // occupation number of an agent state at a particular time
    // time must be between 0 and nTimesteps-1 (final state not implemented at present)
    double operator()(int time, const AGENT &agent) const {
        int beginIndex = Event(time,agent,0);
        assert(beginIndex < size()); // TODO: Implement final state occupation numbers
        int endIndex = beginIndex + AGENT::actDomainSize();
        double occupation = 0.0;
        int eventId;
        for(int eventId=beginIndex; eventId<endIndex; ++eventId) {
            occupation += (*this)[eventId];
        }
        return occupation;
    }

    // time slice for time is in [0...nTimesteps]
    // i.e. final state is implemented
    ModelState<AGENT> operator()(int time) const {
//        glp::SparseVec sparseThis(*this);
//        std::cout << "Getting timeslice from trajectory:" << glp::SparseVec(*this) << std::endl;
        ModelState<AGENT> timeslice;
        int beginIndex = Event(time,AGENT(0),0);
        if(beginIndex >= size()) { // must be final timestep: use act consequence count
            beginIndex = Event(time-1,AGENT(0),0);
            assert(beginIndex < size()); // previous timestep must be in-range
            int endIndex = beginIndex + AGENT::actDomainSize() * AGENT::domainSize();
            for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
                if (double occupation = (*this)[eventId]; fabs(occupation) > tol) {
                    for(const AGENT &consequence : Event<AGENT>(eventId).consequences()) {
                        timeslice[consequence] += occupation;
                    }
                }
            }
        } else {
            int endIndex = beginIndex + AGENT::actDomainSize() * AGENT::domainSize();
            for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
                if (double occupation = (*this)[eventId]; fabs(occupation) > tol)
                    timeslice[Event<AGENT>(eventId).agent()] += occupation;
            }
        }
        return timeslice;
    }

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
        for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
            AGENT agent(agentId);
            std::vector<double> actPMF = agent.timestep(t0State, 0.0);
            assert(occupation <= actPMF.size());
            for(int nthAgent=0; nthAgent < t0State[agentId]; ++nthAgent) {
                typename AGENT::Act act = Random::chooseFromPMF(actPMF);
                trajectory[Event(t,agent,act)] = 1.0;
                t1State += agent.consequences(act);
                actPMF[act] = 0.0;
            }
        }
        t0State.setToZero();
        t0State.swap(t1State);
    }
    return trajectory;
}


#endif //GLPKTEST_TRAJECTORY_H
