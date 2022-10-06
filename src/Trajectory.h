//
// Created by daniel on 04/05/2021.
//

#ifndef GLPKTEST_TRAJECTORY_H
#define GLPKTEST_TRAJECTORY_H

#include <cfloat>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include "ModelState.h"
#include "include/Random.h"
#include "include/StlStream.h"
#include "ABM.h"
#include "EqualityConstraints.h"

template<typename AGENT, int NTIMESTEPS>
class Trajectory: public std::vector<int> {
public:
//    typedef ABM::occupation_type occupation_type;
    typedef AGENT agent_type;

    static constexpr size_t nTimesteps = NTIMESTEPS;
    static constexpr size_t dimension =  AGENT::domainSize*AGENT::actDomainSize*NTIMESTEPS;

    Trajectory(): Trajectory(dimension) { }

//    Trajectory(const ModelState<AGENT> &startState): Trajectory(dimension) {
//        ModelState<AGENT> t0State = startState;
//        ModelState<AGENT> t1State;
//        for (int t = 0; t < nTimesteps; ++t) {
//            for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
//                AGENT agent(agentId);
//                int nAgents = t0State[agentId];
//                for (int actId = 0; actId < AGENT::actDomainSize; ++actId) {
//                    (*this)[indexOf(Event<AGENT>(t, agent, actId))] = 0;
//                }
//                std::vector<double> actPMF = agent.timestep(t0State);
//                for (int a = 0; a < nAgents; ++a) {
//                    int nextAct = Random::nextIntFromDiscrete(actPMF);
//                    (*this)[indexOf(Event<AGENT>(t, agent, nextAct))] += 1;
//                    t1State += agent.consequences(nextAct);
//                }
//            }
//            t0State.setToZero();
//            t0State.swap(t1State);
//        }
//    }


protected:
    Trajectory(size_t nElements): std::vector<value_type>(nElements,0) {}

public:


    const value_type &operator[](const Event<AGENT> &event) const {
        return std::vector<value_type>::operator[](event.id);
    }


    value_type &operator[](int index) {
        return std::vector<value_type>::operator[](index);
    }

    const value_type &operator[](int index) const {
        return std::vector<value_type>::operator[](index);
    }


    static int indexOf(const Event<AGENT> &event) {
        return event.id;
    }

    static EqualityConstraints<value_type> constraints(int nTimesteps) {
        EqualityConstraints<value_type> constraints;
        for(int time = 1; time < nTimesteps; ++time) {
            for(int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
                SparseVec<value_type> coefficients;
                // outgoing edges
                for (int act = 0; act < AGENT::actDomainSize; ++act) {
                    coefficients.insert(indexOf(Event<AGENT>(time, agentState, act)), 1);
                }
                // incoming edges
                for (const Event<AGENT> &inEdge: State<AGENT>::incomingEventsByState[agentState]) {
                    coefficients.insert(indexOf(Event<AGENT>(time-1,inEdge.agent(),inEdge.act())), -1);
                }
                constraints.emplace_back(coefficients,0);
            }
        }
        return constraints;
    }

    friend std::ostream &operator <<(std::ostream &out, const Trajectory<AGENT,NTIMESTEPS> &trajectory) {
        if(trajectory.size() < 2048) { // print whole trajectory
            for (int t = 0; t < NTIMESTEPS; ++t) {
                for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
                    for (int actId = 0; actId < AGENT::actDomainSize; ++actId) {
                        out << trajectory[Event<AGENT>(t, agentId, actId)];
                    }
                    out << " ";
                }
                out << std::endl;
            }
        } else { // print ellipsis
            out << std::vector<Trajectory<AGENT,NTIMESTEPS>::value_type>(trajectory.begin(), trajectory.begin() + 2047) << "..." << std::endl;
        }
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & static_cast<std::vector<value_type> &>(*this);
    }


    //    value_type &operator[](const Event<AGENT> &event) {
//        return std::vector<value_type>::operator[](event.id);
//    }

//    explicit Trajectory(std::vector<value_type> &&rvalue): std::vector<value_type>(rvalue) {
//        assert((this->size()-1)%(AGENT::domainSize*AGENT::actDomainSize) == 0);
//    }
//
//    explicit Trajectory(const std::vector<value_type> &lvalue, int nTimesteps): std::vector<value_type>(lvalue) {
//        resize(nTimesteps*AGENT::domainSize*AGENT::actDomainSize,0);
//    }

//    ModelState<AGENT> modelState(int time) const { return ModelState<AGENT>(*this, time); }

    // occupation number of an agent state at a particular time
//    value_type operator[](const State<AGENT> &state) const {
//        assert(state.time >= 0 && state.time <= nTimesteps());
//        return(state.time < nTimesteps())?state.forwardOccupation(*this):state.backwardOccupation(*this);
//    }



//    static int nTimesteps() { return NTIMESTEPS; }  //{  return size()/(AGENT::domainSize*AGENT::actDomainSize); }

//    int dimension() const { return size(); }
//
//    static int dimension(int nTimesteps) { return AGENT::domainSize*AGENT::actDomainSize*nTimesteps; }

//    const ModelState<AGENT> &temporaryPartialModelState(int time, const std::vector<int> &agentIds) const {
//        static thread_local ModelState<AGENT> state;
//        for(int agentId: agentIds)
//            state[agentId] = (*this)[State<AGENT>(time, agentId)];
//        return state;
//    }

//    const ModelState<AGENT> &temporaryPartialModelState(int time, const std::vector<AGENT> &agentIds) const {
//        static thread_local ModelState<AGENT> state;
//        for(AGENT agent: agentIds)
//            state[agent] = (*this)[State<AGENT>(time, agent)];
//        return state;
//    }


//    Trajectory<AGENT> slice(int fromTimestep, int nTimesteps) const {
//        Trajectory<AGENT> slice(nTimesteps);
//        int beginIndex = Event(fromTimestep,AGENT(0),0);
//        int endIndex = Event(fromTimestep+nTimesteps,AGENT(0),0);
//        int sliceBeginIndex = Event(0,AGENT(0),0);
//        std::copy(begin()+beginIndex, begin()+endIndex, slice.begin()+sliceBeginIndex);
//        return slice;
//    }


    // gives the coefficients to create the givven state
//    static SparseVec<value_type> coefficients(const State<AGENT> &state) {
//        SparseVec<value_type> coeffs;
//        coeffs.indices = state.forwardOccupationDependencies();
//        coeffs.values.resize(coeffs.indices.size(),1);
//        return coeffs;
//    }

};


#endif //GLPKTEST_TRAJECTORY_H
