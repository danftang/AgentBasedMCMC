//
// Created by daniel on 11/05/2021.
//

#ifndef GLPKTEST_TIMESTEP_H
#define GLPKTEST_TIMESTEP_H

#include "glpkpp.h"

template<typename AGENT>
class Timestep: public glp::SparseVec {
public:
    Timestep() { }
    Timestep(glp::SparseVec &&rvalue): glp::SparseVec(0) {
        swap(rvalue);
    }

    Timestep(const glp::SparseVec &lvalue): glp::SparseVec(lvalue) { }

    Timestep &operator =(glp::SparseVec &&rvalue) {
        swap(rvalue);
        return *this;
    }

    Timestep() &operator =(const glp::SparseVec &lvalue) {
        glp::SparseVec::operator=(lvalue);
        return *this;
    }

    // event count
    double operator()(const AGENT &agent, const typename AGENT::Act &act) const {
        return (*this)[Event(0,agent,act)];
    }

    // occupation number
    double operator[](const AGENT &agent) const {
        int beginIndex = Event(0,agent,typename AGENT::Act(0));
        int endIndex = beginIndex + AGENT::actDomainSize();
        double occupation = 0.0;
        int index;
        for(int i=1; i<=sparseSize(); ++i) {
            index = indices[i];
            if(beginIndex <= index && index < endIndex) occupation += values[i];
        }
        return occupation;
    }

    friend std::ostream &operator <<(std::ostream &out, const Timestep<AGENT> &timestep) {
        for(auto eventCount: timestep) {
            Event<AGENT> event = eventCount.index;
            if(event.time() != timestep) {
                out << "time = " << event.time() << std::endl;
                timestep = event.time();
            }
            out << "  " << eventCount.value << "*" << event << std::endl;
        }

        return out;
    }
};


#endif //GLPKTEST_TIMESTEP_H
