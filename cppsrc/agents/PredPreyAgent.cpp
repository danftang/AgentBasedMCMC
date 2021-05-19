//
// Created by daniel on 18/05/2021.
//

#include <cassert>
#include "PredPreyAgent.h"
#include "../State.h"

int PredPreyAgent::GRIDSIZE = 16;

std::vector<PredPreyAgent> PredPreyAgent::consequences(PredPreyAgent::Act act) const {
    switch(act) {
        case DIE:
            return std::vector<PredPreyAgent>();
        case MOVELEFT:
            return std::vector<PredPreyAgent>( { PredPreyAgent(xLeft(),yPosition(),type()) } );
        case MOVERIGHT:
            return std::vector<PredPreyAgent>( { PredPreyAgent(xRight(),yPosition(),type()) } );
        case MOVEUP:
            return std::vector<PredPreyAgent>( { PredPreyAgent(xPosition(),yUp(),type()) } );
        case MOVEDOWN:
            return std::vector<PredPreyAgent>( { PredPreyAgent(xPosition(),yDown(),type()) } );
        case GIVEBIRTH:
            return std::vector<PredPreyAgent>( { *this, PredPreyAgent(xRight(),yPosition(),type()) } );
        default:
            assert(false);
    }
}

std::vector<double> PredPreyAgent::timestep(const ModelState<PredPreyAgent> &others) const {
        static double pPredBirthGivenPrey = 0.5; // birth prob given prey
        static double pPredDie = 0.07; // death prob
        static double pPreyBirth = 0.06; // birth prob
        static double pPreyDie = 0.03; // death prob
        static double pPreyEatenGivenPred = 0.55; // death prob given pred

        std::vector<double> actDistribution(actDomainSize());
        if(type() == PREDATOR) {
            actDistribution[DIE] = pPredDie;
            if(surroundingCountOf(PREY, others) >= 1.0) {
                actDistribution[GIVEBIRTH] = pPredBirthGivenPrey;
            }
        } else {
            actDistribution[GIVEBIRTH] = pPreyBirth;
            actDistribution[DIE] = pPreyDie;
            if(surroundingCountOf(PREDATOR, others) >= 1.0) {
                actDistribution[DIE] += pPreyEatenGivenPred;
            }
        }
        double moveProb = 0.25*(1.0 - actDistribution[DIE] - actDistribution[GIVEBIRTH]);
        actDistribution[MOVEUP] = moveProb;
        actDistribution[MOVEDOWN] = moveProb;
        actDistribution[MOVELEFT] = moveProb;
        actDistribution[MOVERIGHT] = moveProb;
        return actDistribution;
}

std::vector<glp::Constraint> PredPreyAgent::constraints(int time, PredPreyAgent::Act act) const {
    if(type() == PREDATOR && act == GIVEBIRTH) {
        return std::vector<glp::Constraint>({
            1.0*State(time,PredPreyAgent(xLeft(),yPosition(), PREY)) +
            1.0*State(time,PredPreyAgent(xRight(),yPosition(), PREY)) +
            1.0*State(time,PredPreyAgent(xPosition(),yUp(),PREY)) +
            1.0*State(time,PredPreyAgent(xPosition(),yDown(),PREY)) >= 1.0
        });
    }
    return std::vector<glp::Constraint>();
}
