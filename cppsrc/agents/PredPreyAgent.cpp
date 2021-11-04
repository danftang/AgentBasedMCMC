//
// Created by daniel on 18/05/2021.
//

//#include <cassert>
//#include "PredPreyAgent.h"
//#include "../State.h"
//
//
//template<int GRIDSIZE>
//std::vector<PredPreyAgent<GRIDSIZE>> PredPreyAgent<GRIDSIZE>::consequences(PredPreyAgent<GRIDSIZE>::Act act) const {
//    switch(act) {
//        case DIE:
//            return std::vector<PredPreyAgent<GRIDSIZE>>();
//        case MOVELEFT:
//            return std::vector<PredPreyAgent<GRIDSIZE>>( { PredPreyAgent<GRIDSIZE>(xLeft(),yPosition(),type()) } );
//        case MOVERIGHT:
//            return std::vector<PredPreyAgent<GRIDSIZE>>( { PredPreyAgent<GRIDSIZE>(xRight(),yPosition(),type()) } );
//        case MOVEUP:
//            return std::vector<PredPreyAgent<GRIDSIZE>>( { PredPreyAgent<GRIDSIZE>(xPosition(),yUp(),type()) } );
//        case MOVEDOWN:
//            return std::vector<PredPreyAgent<GRIDSIZE>>( { PredPreyAgent<GRIDSIZE>(xPosition(),yDown(),type()) } );
//        case GIVEBIRTH:
//            return std::vector<PredPreyAgent<GRIDSIZE>>( { *this, PredPreyAgent<GRIDSIZE>(xRight(),yPosition(),type()) } );
//    }
//    assert(false); // unrecognized act
//    return std::vector<PredPreyAgent<GRIDSIZE>>();
//}
//
//
//// If allowInfeasibleActs is true, then if an act is infeasible then its probability is not set to zero but
//// is the expectation averaged over the trajectoryPrior distribution of 'others' (note that this has the
//// consequence that the sum of the acts no longer adds to 1 since different acts are relative to different
//// distributions of 'others'). This is to enable different uses: for a forward run we set allowInfeasibleActs to
//// false and get a LogPMF over acts, whereas when calculating the probability of a exactEndState by setting allowInfeasibleActs
//// we extend probability calculation to infeasible trajectories in a reasonably smooth manner.
//template<int GRIDSIZE>
//std::vector<double> PredPreyAgent<GRIDSIZE>::timestep(const ModelState<PredPreyAgent<GRIDSIZE>> &others) const {
//
//    std::vector<double> actDistribution(actDomainSize(),0.0);
//    if (type() == PREDATOR) {
//        actDistribution[DIE] = pPredDie;
//        if (surroundingCountOf(PREY, others) >= 1.0) {
//            actDistribution[GIVEBIRTH] = pPredBirthGivenPrey;
//        }
//    } else {
//        actDistribution[GIVEBIRTH] = pPreyBirth;
//        actDistribution[DIE] = pPreyDie;
//        if (surroundingCountOf(PREDATOR, others) >= 1.0) {
//            actDistribution[DIE] += pPreyEatenGivenPred;
//        }
//    }
//    double moveProb = 0.25 * (1.0 - actDistribution[DIE] - actDistribution[GIVEBIRTH]);
//    actDistribution[MOVEUP] = moveProb;
//    actDistribution[MOVEDOWN] = moveProb;
//    actDistribution[MOVELEFT] = moveProb;
//    actDistribution[MOVERIGHT] = moveProb;
//
//    return actDistribution;
//}
//
//
//template<int GRIDSIZE>
//std::vector<double> PredPreyAgent<GRIDSIZE>::marginalTimestep() const {
//    constexpr double pPreyNeighbourGivenPred = 0.2;
//    constexpr double pPredNeighbourGivenPrey = 0.1;
//
//    std::vector<double> actDistribution(actDomainSize(),0.0);
//    if (type() == PREDATOR) {
//        actDistribution[DIE] = pPredDie;
//        actDistribution[GIVEBIRTH] = pPredBirthGivenPrey*pPreyNeighbourGivenPred;
//    } else {
//        actDistribution[GIVEBIRTH] = pPreyBirth;
//        actDistribution[DIE] = pPreyDie + pPreyEatenGivenPred*pPredNeighbourGivenPrey;
//    }
//    double moveProb = 0.25 * (1.0 - actDistribution[DIE] - actDistribution[GIVEBIRTH]);
//    actDistribution[MOVEUP] = moveProb;
//    actDistribution[MOVEDOWN] = moveProb;
//    actDistribution[MOVELEFT] = moveProb;
//    actDistribution[MOVERIGHT] = moveProb;
//
//    return actDistribution;
//}
//
//
//template<int GRIDSIZE>
//std::vector<glp::Constraint> PredPreyAgent<GRIDSIZE>::constraints(int time, PredPreyAgent<GRIDSIZE>::Act act) const {
//    if(type() == PREDATOR && act == GIVEBIRTH) {
//        return std::vector<glp::Constraint>({
//            1.0*State(time,PredPreyAgent<GRIDSIZE>(xLeft(),yPosition(), PREY)) +
//            1.0*State(time,PredPreyAgent<GRIDSIZE>(xRight(),yPosition(), PREY)) +
//            1.0*State(time,PredPreyAgent<GRIDSIZE>(xPosition(),yUp(),PREY)) +
//            1.0*State(time,PredPreyAgent<GRIDSIZE>(xPosition(),yDown(),PREY)) >= 1.0
//        });
//    }
//    return std::vector<glp::Constraint>();
//}
