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

    // execute forward from start state
    Trajectory(int nTimesteps, const ModelState<AGENT> &startState): Trajectory(nTimesteps) {
        ModelState<AGENT> t0State = startState;
        ModelState<AGENT> t1State;
        for(int t=0; t<nTimesteps; ++t) {
            for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
                AGENT agent(agentId);
                std::vector<double> actPMF = agent.timestep(t0State, 0.0);
                assert(t0State[agentId] <= actPMF.size());
                for(int nthAgent=0; nthAgent < t0State[agentId]; ++nthAgent) {
                    typename AGENT::Act act = Random::chooseFromPMF(actPMF);
                    operator[](Event(t,agent,act)) = 1.0;
                    t1State += agent.consequences(act);
                    actPMF[act] = 0.0;
                }
            }
            t0State.setToZero();
            t0State.swap(t1State);
        }
    }


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
        assert(time>=0 && time<=nTimesteps());
        if(time == nTimesteps()) return endState();
        ModelState<AGENT> timeslice;
        int beginIndex = Event(time,AGENT(0),0);
        int endIndex = beginIndex + AGENT::actDomainSize() * AGENT::domainSize();
        for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
            if (double occupation = (*this)[eventId]; fabs(occupation) > tol)
                timeslice[Event<AGENT>(eventId).agent()] += occupation;
        }
        return timeslice;
    }


    ModelState<AGENT> endState() const {
        ModelState<AGENT> timeslice;
        int beginIndex = Event(nTimesteps(),AGENT(0),0);
        int endIndex = beginIndex + AGENT::actDomainSize() * AGENT::domainSize();
        for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
            if (double occupation = (*this)[eventId]; fabs(occupation) > tol) {
                for(const AGENT &consequence : Event<AGENT>(eventId).consequences()) {
                    timeslice[consequence] += occupation;
                }
            }
        }
        return timeslice;
    }


    int nTimesteps() const { return (size()-1)/(AGENT::domainSize()*AGENT::actDomainSize()); }



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

#endif //GLPKTEST_TRAJECTORY_H
