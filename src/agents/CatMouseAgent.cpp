//
// Created by daniel on 26/04/2021.
//

#include "CatMouseAgent.h"

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


double CatMouseAgent::marginalTimestep(Act act) const {
    if (type() == CAT) {
        if(act == MOVE) return pCatMove; else return 1.0 - pCatMove;
    }
    return 1.0;
}


std::vector<CatMouseAgent> CatMouseAgent::consequences(Act act) const {
    if(act == MOVE) {
        return std::vector<CatMouseAgent>({ CatMouseAgent(type(), Position((position()+1)%2)) });
    } else {
        return std::vector<CatMouseAgent>({ CatMouseAgent(type(), position()) });
    }
}


// result of static analysis of timestep member function...
std::vector<Constraint<ABM::occupation_type>> CatMouseAgent::constraints(int time, Act act) const {
    if(type() == MOUSE) {
        if(act == MOVE) {
            return std::vector({ 1*State(time,CatMouseAgent(CAT, position())) >= 1 });
        } else {
            return std::vector({ 1*State(time,CatMouseAgent(CAT, position())) <= 0 });
        }
    } else {
        return std::vector<Constraint<ABM::occupation_type>>();
    }
}

std::vector<CatMouseAgent> CatMouseAgent::neighbours() {
    if(type() == CAT) return { CatMouseAgent(MOUSE, position()) };
    return {};
}
