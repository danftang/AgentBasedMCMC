//
// Created by daniel on 04/05/2021.
//

#ifndef GLPKTEST_TRAJECTORY_H
#define GLPKTEST_TRAJECTORY_H

#include <cfloat>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include "ModelState.h"
#include "StateTrajectory.h"
#include "include/Random.h"
#include "ActFermionicDistribution.h"
#include "ABM.h"
#include "EqualityConstraints.h"

template<typename AGENT>
class Trajectory: public std::vector<ABM::occupation_type> {
public:
    typedef ABM::occupation_type occupation_type;
    typedef AGENT agent_type;

    Trajectory() = default;

    explicit Trajectory(int nTimesteps): std::vector<value_type>(dimension(nTimesteps)) { }

    explicit Trajectory(std::vector<value_type> &&rvalue): std::vector<value_type>(rvalue) {
        assert((this->size()-1)%(AGENT::domainSize()*AGENT::actDomainSize()) == 0);
    }

    explicit Trajectory(const std::vector<value_type> &lvalue, int nTimesteps): std::vector<value_type>(lvalue) {
        resize(nTimesteps*AGENT::domainSize()*AGENT::actDomainSize(),0);
    }

    ModelState<AGENT> modelState(int time) const { return ModelState<AGENT>(*this, time); }

    // occupation number of an agent state at a particular time
    value_type operator[](const State<AGENT> &state) const {
        assert(state.time >= 0 && state.time <= nTimesteps());
        return(state.time < nTimesteps())?state.forwardOccupation(*this):state.backwardOccupation(*this);
    }

    const value_type &operator[](const Event<AGENT> &event) const {
        return std::vector<value_type>::operator[](event.id);
    }

    value_type &operator[](const Event<AGENT> &event) {
        return std::vector<value_type>::operator[](event.id);
    }

    value_type &operator[](int index) {
        return std::vector<value_type>::operator[](index);
    }

    const value_type &operator[](int index) const {
        return std::vector<value_type>::operator[](index);
    }



    int nTimesteps() const { return size()/(AGENT::domainSize()*AGENT::actDomainSize()); }

    int dimension() const { return size(); }

    static int dimension(int nTimesteps) { return AGENT::domainSize()*AGENT::actDomainSize()*nTimesteps; }

//    const ModelState<AGENT> &temporaryPartialModelState(int time, const std::vector<int> &agentIds) const {
//        static thread_local ModelState<AGENT> state;
//        for(int agentId: agentIds)
//            state[agentId] = (*this)[State<AGENT>(time, agentId)];
//        return state;
//    }

    const ModelState<AGENT> &temporaryPartialModelState(int time, const std::vector<AGENT> &agentIds) const {
        static thread_local ModelState<AGENT> state;
        for(AGENT agent: agentIds)
            state[agent] = (*this)[State<AGENT>(time, agent)];
        return state;
    }


    Trajectory<AGENT> slice(int fromTimestep, int nTimesteps) const {
        Trajectory<AGENT> slice(nTimesteps);
        int beginIndex = Event(fromTimestep,AGENT(0),0);
        int endIndex = Event(fromTimestep+nTimesteps,AGENT(0),0);
        int sliceBeginIndex = Event(0,AGENT(0),0);
        std::copy(begin()+beginIndex, begin()+endIndex, slice.begin()+sliceBeginIndex);
        return slice;
    }


    // gives the coefficients to create the givven state
    static SparseVec<value_type> coefficients(const State<AGENT> &state) {
        SparseVec<value_type> coeffs;
        coeffs.indices = state.forwardOccupationDependencies();
        coeffs.values.resize(coeffs.indices.size(),1);
        return coeffs;
    }

    static int indexOf(const Event<AGENT> &event) {
        return event.id;
    }

    static EqualityConstraints<value_type> constraints(int nTimesteps) {
        EqualityConstraints<value_type> constraints;
        for(int time = 1; time < nTimesteps; ++time) {
            for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                SparseVec<ABM::coefficient_type> coefficients;
                // outgoing edges
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
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



private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & static_cast<std::vector<value_type> &>(*this);
    }
};


#endif //GLPKTEST_TRAJECTORY_H
