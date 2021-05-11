//
// Created by daniel on 26/04/2021.
//

#include "CatMouseAgent.h"
#include "State.h"

// returns PMF over acts
std::vector<double> CatMouseAgent::timestep(std::map<CatMouseAgent,double> others) {
    static constexpr double pCatMove = 0.25;
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

std::vector<CatMouseAgent> CatMouseAgent::consequences(Act act) {
    if(act == MOVE) {
        return std::vector<CatMouseAgent>({ CatMouseAgent(type(), Position((position()+1)%2)) });
    } else {
        return std::vector<CatMouseAgent>({ CatMouseAgent(type(), position()) });
    }
}


// result of static analysis of timestep member function...
std::vector<glp::Constraint> CatMouseAgent::constraints(int time, Act act) {
    if(type() == MOUSE) {
        if(act == MOVE) {
            return std::vector({ 1.0*State(time,CatMouseAgent(CAT, position())) >= 1 });
        } else {
            return std::vector({ 1.0*State(time,CatMouseAgent(CAT, position())) == 0 });
        }
    } else {
        return std::vector<glp::Constraint>();
    }
}
