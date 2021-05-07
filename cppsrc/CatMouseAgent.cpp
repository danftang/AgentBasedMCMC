//
// Created by daniel on 26/04/2021.
//

#include "CatMouseAgent.h"
#include "State.h"

// returns PMF over acts
std::vector<double> CatMouseAgent::timestep(std::multiset<CatMouseAgent> others) {

//        val pCatMove = 0.25
//
//        return if(type == AgentType.CAT) {
//            arrayOf(pCatMove, 1.0-pCatMove)
//        } else {
//            if(others[CatMouseAgent(AgentType.CAT, position)] >= 1) {
//                arrayOf(1.0,0.0)
//            } else {
//                arrayOf(0.0,1.0)
//            }
//        }
    return std::vector<double>(); // placeholder
}

std::vector<CatMouseAgent> CatMouseAgent::consequences(Act act) {
    if(act == MOVE) {
        return std::vector<CatMouseAgent>({ CatMouseAgent(type(), Position((position()+1)%2)) });
    } else {
        return std::vector<CatMouseAgent>({ CatMouseAgent(type(), position()) });
    }
}

std::vector<glp::Constraint> CatMouseAgent::constraints(Act act) {
    if(type() == MOUSE) {
        if(act == MOVE) {
            return std::vector({ 1.0*State(0,CatMouseAgent(CAT, position())) >= 1 });
        } else {
            return std::vector({ 1.0*State(0,CatMouseAgent(CAT, position())) == 0 });
        }
    } else {
        return std::vector<glp::Constraint>();
    }
}
