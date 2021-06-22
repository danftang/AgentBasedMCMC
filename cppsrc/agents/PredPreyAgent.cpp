//
// Created by daniel on 18/05/2021.
//

#include <cassert>
#include "PredPreyAgent.h"
#include "../State.h"

int PredPreyAgent::GRIDSIZE = 16; // default gridsize


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


// If allowInfeasibleActs is true, then if an act is infeasible then its probability is not set to zero but
// is the expectation averaged over the prior distribution of 'others' (note that this has the
// consequence that the sum of the acts no longer adds to 1 since different acts are relative to different
// distributions of 'others'). This is to enable different uses: for a forward run we set allowInfeasibleActs to
// false and get a PMF over acts, whereas when calculating the probability of a trajectory by setting allowInfeasibleActs
// we extend probability calculation to infeasible trajectories in a reasonably smooth manner.
std::vector<double> PredPreyAgent::timestep(const ModelState<PredPreyAgent> &others, double infeasibilityProb) const {

    std::vector<double> actDistribution(actDomainSize(),0.0);
    if (type() == PREDATOR) {
        actDistribution[DIE] = pPredDie;
        if (surroundingCountOf(PREY, others) >= 1.0) {
            actDistribution[GIVEBIRTH] = pPredBirthGivenPrey;
        }
    } else {
        actDistribution[GIVEBIRTH] = pPreyBirth;
        actDistribution[DIE] = pPreyDie;
        if (surroundingCountOf(PREDATOR, others) >= 1.0) {
            actDistribution[DIE] += pPreyEatenGivenPred;
        }
    }
    double moveProb = 0.25 * (1.0 - actDistribution[DIE] - actDistribution[GIVEBIRTH]);
    actDistribution[MOVEUP] = moveProb;
    actDistribution[MOVEDOWN] = moveProb;
    actDistribution[MOVELEFT] = moveProb;
    actDistribution[MOVERIGHT] = moveProb;

    if(actDistribution[GIVEBIRTH] == 0.0) actDistribution[GIVEBIRTH] = infeasibilityProb;
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
