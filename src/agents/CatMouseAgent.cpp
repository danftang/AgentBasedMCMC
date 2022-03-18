//
// Created by daniel on 26/04/2021.
//

#include "CatMouseAgent.h"
#include "../State.h"

// returns LogPMF over acts
std::vector<double> CatMouseAgent::timestep(const ModelState<CatMouseAgent> &others) const {
    std::vector<double> actPmf(actDomainSize());

    if (type() == CAT) {
        actPmf[MOVE] = pCatMove;
        actPmf[STAYPUT] = 1.0 - pCatMove;
    } else {
        if (others[CatMouseAgent(CAT, position())] >= 1.0) {
            actPmf[MOVE] = 1.0;
            actPmf[STAYPUT] = 0.0;
        } else {
            actPmf[MOVE] = 0.0;
            actPmf[STAYPUT] = 1.0;
        }
    }
    return actPmf;
}


std::vector<double> CatMouseAgent::marginalTimestep() const {
    std::vector<double> actPmf(actDomainSize());

    if (type() == CAT) {
        actPmf[MOVE] = pCatMove;
        actPmf[STAYPUT] = 1.0 - pCatMove;
    } else {
        actPmf[MOVE] = 1.0;     // given that all constraints are satisfied
        actPmf[STAYPUT] = 1.0;
    }
    return actPmf;
}


std::vector<CatMouseAgent> CatMouseAgent::consequences(Act act) const {
    if(act == MOVE) {
        return std::vector<CatMouseAgent>({ CatMouseAgent(type(), Position((position()+1)%2)) });
    } else {
        return std::vector<CatMouseAgent>({ CatMouseAgent(type(), position()) });
    }
}


// result of static analysis of timestep member function...
std::vector<glp::Constraint> CatMouseAgent::constraints(int time, Act act) const {
    if(type() == MOUSE) {
        if(act == MOVE) {
            return std::vector({ 1.0*State(time,CatMouseAgent(CAT, position())) >= 1 });
        } else {
            return std::vector({ 1.0*State(time,CatMouseAgent(CAT, position())) <= 0 });
        }
    } else {
        return std::vector<glp::Constraint>();
    }
}
