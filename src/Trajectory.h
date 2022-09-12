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
#include "Random.h"
#include "ActFermionicDistribution.h"
#include "ABM.h"

template<typename AGENT>
class Trajectory: public std::vector<ABM::occupation_type> {
public:
    typedef ABM::occupation_type occupation_type;
    typedef AGENT agent_type;


    explicit Trajectory(int nTimesteps): std::vector<value_type>(dimension(nTimesteps)) { }

    explicit Trajectory(std::vector<value_type> &&rvalue): std::vector<value_type>(rvalue) {
        assert((this->size()-1)%(AGENT::domainSize()*AGENT::actDomainSize()) == 0);
    }

    explicit Trajectory(const std::vector<value_type> &lvalue, int nTimesteps): std::vector<value_type>(lvalue) {
        resize(nTimesteps*AGENT::domainSize()*AGENT::actDomainSize(),0);
    }


    // occupation number of an agent state at a particular time
    value_type operator[](const State<AGENT> &state) const {
        assert(state.time >= 0 && state.time <= nTimesteps());
        return(state.time < nTimesteps())?state.forwardOccupation(*this):state.backwardOccupation(*this);
    }

    value_type &operator[](const Event<AGENT> &event) {
        return std::vector<value_type>::operator[](event.id);
    }

    const value_type &operator[](const Event<AGENT> &event) const {
        return std::vector<value_type>::operator[](event.id);
    }

    value_type &operator[](int eventId) {
        return std::vector<value_type>::operator[](eventId);
    }

    const value_type &operator[](int eventId) const {
        return std::vector<value_type>::operator[](eventId);
    }


    ModelState<AGENT> endState() const {
        return ModelState<AGENT>(*this, nTimesteps(), nTimesteps());
    }

    int nTimesteps() const { return size()/(AGENT::domainSize()*AGENT::actDomainSize()); }

    int dimension() const { return size(); }

    static int dimension(int nTimesteps) { return AGENT::domainSize()*AGENT::actDomainSize()*nTimesteps; }


    Trajectory<AGENT> slice(int fromTimestep, int nTimesteps) const {
        Trajectory<AGENT> slice(nTimesteps);
        int beginIndex = Event(fromTimestep,AGENT(0),0);
        int endIndex = Event(fromTimestep+nTimesteps,AGENT(0),0);
        int sliceBeginIndex = Event(0,AGENT(0),0);
        std::copy(begin()+beginIndex, begin()+endIndex, slice.begin()+sliceBeginIndex);
        return slice;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & static_cast<std::vector<value_type> &>(*this);
    }
};


#endif //GLPKTEST_TRAJECTORY_H
