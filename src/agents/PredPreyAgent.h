//
// Created by daniel on 18/05/2021.
//

#ifndef GLPKTEST_PREDPREYAGENT_H
#define GLPKTEST_PREDPREYAGENT_H

#include <vector>
#include <map>
#include "../ModelState.h"
#include "../Constraint.h"

class PredPreyAgentBase {
public:
    typedef int Act;

    enum ActNames {
        GIVEBIRTH,
        DIE,
        MOVELEFT,
        MOVERIGHT,
        MOVEUP,
        MOVEDOWN
    };


    enum Type {
        PREDATOR,
        PREY
    };


    static constexpr double pPred = 0.05; // steady state probability that there are no predators on a square
    static constexpr double pPrey = 0.05; // steady state probability that there are no prey on a square
    static constexpr double pNoPred = 1.0-pPred; // steady state probability that there are no predators on a square
    static constexpr double pNoPrey = 1.0-pPrey; // steady state probability that there are no prey on a square
    static constexpr double clustering = 0.33; // tendency to cluster / instability
    static constexpr double pNoPredClose = pNoPred*pNoPred*pNoPred*pNoPred; // no pow function for constexpr :-(
    static constexpr double pNoPreyClose = pNoPrey*pNoPrey*pNoPrey*pNoPrey;
    static constexpr double pPredClose = 1.0 - pNoPredClose; // no pow function for constexpr :-(
    static constexpr double pPreyClose = 1.0 - pNoPreyClose;
    // values for which there is a steady state prior of uniform, uncorrelated Poisson distribution
    // At steady state probability of birth must equal prob of death, so:
    //      pPredDie = pPreyClose*pPredBirthGivenPrey
    //      pPreyBirth = pPreyDie + pPredClose*pPreyEatenGivenPred
//    static constexpr double pPredBirthGivenPrey = clustering;
//    static constexpr double pPredDie = pPreyClose*pPredBirthGivenPrey;
//    static constexpr double pPreyEatenGivenPred = clustering;
//    static constexpr double pPreyDie = 0.1; // stability/mixing?
//    static constexpr double pPreyBirth = pPreyDie + pPredClose*pPreyEatenGivenPred;

    // If we fix the probability of movement, so that the pBirth + pDeath is independent of surrounding
    // and assume no pred death in the presence of prey then we have
    // P(predDie|noPreyClose) = c
    // P(predDie|preyClose) = P(predBorn|preyClose)P(preyClose)/c - P(noPreyClose)
    // P(predBorn|preyClose) = c - P(predDie|preyClose)
    //
    // And if we assume no prey birth in the presence of preds then
    // P(predDie|predClose)P(predClose) = (cP(predClose)  + cP(noPredClose))/2
    // P(predBirth|predClose) = c - P(predDie|predClose)
    static const double lpPredBirthGivenPrey;
    static const double lpPredDeathGivenPrey;
    static const double lpPredBirthGivenNoPrey;
    static const double lpPredDeathGivenNoPrey;

    static const double lpPreyBirthGivenPred;
    static const double lpPreyDeathGivenPred;
    static const double lpPreyBirthGivenNoPred;
    static const double lpPreyDeathGivenNoPred;

    static const double lpMove;

    // Agent Domain stuff
    static constexpr int actDomainSize() { return 6; }
};

inline const double PredPreyAgentBase::lpPredBirthGivenPrey = log(clustering);
inline const double PredPreyAgentBase::lpPredDeathGivenPrey = log(0.0);
inline const double PredPreyAgentBase::lpPredBirthGivenNoPrey = log(0.5 * clustering * (1.0 - pPreyClose / pNoPreyClose));
inline const double PredPreyAgentBase::lpPredDeathGivenNoPrey = log(clustering - exp(lpPredBirthGivenNoPrey));

inline const double PredPreyAgentBase::lpPreyBirthGivenPred = log(0.0);
inline const double PredPreyAgentBase::lpPreyDeathGivenPred = log(clustering);
inline const double PredPreyAgentBase::lpPreyBirthGivenNoPred = log(0.5 * clustering * (1.0 + pPredClose / pNoPredClose));
inline const double PredPreyAgentBase::lpPreyDeathGivenNoPred = log(clustering - exp(lpPreyBirthGivenNoPred));

inline const double PredPreyAgentBase::lpMove = log(0.25 * (1.0 - clustering));


template<int GRIDSIZE>
class PredPreyAgent: public PredPreyAgentBase {
public:
    // Agent Domain stuff
    static constexpr int domainSize() { return GRIDSIZE*GRIDSIZE*2; }

    int stateId;

    PredPreyAgent(int ordinal): stateId(ordinal) {}
    PredPreyAgent( int x, int y, Type type): stateId(x + GRIDSIZE*y + GRIDSIZE*GRIDSIZE*type) { }

//    std::vector<double> timestep(const ModelState<PredPreyAgent<GRIDSIZE>> &others) const;
    std::vector<PredPreyAgent<GRIDSIZE>> consequences(Act act) const; // the consequences of an act

    template<class TRAJECTORY>
    static double logEventProb(const Event<PredPreyAgent<GRIDSIZE>> &event, const TRAJECTORY &trajectory);

    static TrajectoryDependencies<PredPreyAgent<GRIDSIZE>> eventProbDependencies(const Event<PredPreyAgent<GRIDSIZE>> &event);


    operator int() const { return stateId; }
    int xPosition() const { return stateId%GRIDSIZE; }
    int yPosition() const { return (stateId/GRIDSIZE)%GRIDSIZE; }
    Type type() const { return Type(stateId/(GRIDSIZE*GRIDSIZE)); }

    int xLeft() const { return((stateId+GRIDSIZE-1)%GRIDSIZE); }
    int xRight() const { return((stateId+1)%GRIDSIZE); }
    int yUp() const { return((stateId/GRIDSIZE+1)%GRIDSIZE); }
    int yDown() const { return((stateId/GRIDSIZE + GRIDSIZE-1)%GRIDSIZE); }


//    double marginalTimestep(Act act) const; // returns the probability of this act, given an agent in this state
    // returns the constraints implied by the given act of this agent
//    std::vector<Constraint<ABM::occupation_type>> constraints(int time, Act act) const; // to be generated automatically by static analysis...eventually.

    friend std::ostream &operator <<(std::ostream &out, const PredPreyAgent<GRIDSIZE> &agent) {

        out << (agent.type()==PREY?"PREY":"PRED") << ":(" << agent.xPosition() << "," << agent.yPosition() <<")";
        return out;
    }


    // Count of agents of given type north, south, east and west of this agent.
    double surroundingCountOf(Type type, const ModelState<PredPreyAgent<GRIDSIZE>> &others) const {
        return others[PredPreyAgent(xRight(),yPosition(),type)] +
        others[PredPreyAgent(xLeft(),yPosition(),type)] +
        others[PredPreyAgent(xPosition(),yUp(),type)] +
        others[PredPreyAgent(xPosition(),yDown(),type)];
    }

    template<class TRAJECTORY>
    bool hasAnySurrounding(Type type, int time, const TRAJECTORY &trajectory) const {
        return
        trajectory[State(time,PredPreyAgent(xRight(),yPosition(),type))] ||
        trajectory[State(time,PredPreyAgent(xLeft(),yPosition(),type))] ||
        trajectory[State(time,PredPreyAgent(xPosition(),yUp(),type))] ||
        trajectory[State(time,PredPreyAgent(xPosition(),yDown(),type))];
    }



};


template<int GRIDSIZE>
std::vector<PredPreyAgent<GRIDSIZE>> PredPreyAgent<GRIDSIZE>::consequences(PredPreyAgent<GRIDSIZE>::Act act) const {
    switch(act) {
        case DIE:
            return std::vector<PredPreyAgent<GRIDSIZE>>();
        case MOVELEFT:
            return std::vector<PredPreyAgent<GRIDSIZE>>( { PredPreyAgent<GRIDSIZE>(xLeft(),yPosition(),type()) } );
        case MOVERIGHT:
            return std::vector<PredPreyAgent<GRIDSIZE>>( { PredPreyAgent<GRIDSIZE>(xRight(),yPosition(),type()) } );
        case MOVEUP:
            return std::vector<PredPreyAgent<GRIDSIZE>>( { PredPreyAgent<GRIDSIZE>(xPosition(),yUp(),type()) } );
        case MOVEDOWN:
            return std::vector<PredPreyAgent<GRIDSIZE>>( { PredPreyAgent<GRIDSIZE>(xPosition(),yDown(),type()) } );
        case GIVEBIRTH:
            return std::vector<PredPreyAgent<GRIDSIZE>>( { *this, PredPreyAgent<GRIDSIZE>(xRight(),yPosition(),type()) } );
    }
    assert(false); // unrecognized act
    return std::vector<PredPreyAgent<GRIDSIZE>>();
}



template<int GRIDSIZE>
template<class TRAJECTORY>
double PredPreyAgent<GRIDSIZE>::logEventProb(const Event<PredPreyAgent<GRIDSIZE>> &event, const TRAJECTORY &trajectory) {
    switch(event.act()) {
        case MOVEUP:
        case MOVEDOWN:
        case MOVELEFT:
        case MOVERIGHT:
            return lpMove;
        case GIVEBIRTH:
            if (event.agent().type() == PREDATOR) {
                return (event.agent().hasAnySurrounding(PREY, event.time(), trajectory)) ? lpPredBirthGivenPrey : lpPredBirthGivenNoPrey;
            } else {
                return (event.agent().hasAnySurrounding(PREDATOR, event.time(), trajectory)) ? lpPreyBirthGivenPred : lpPreyBirthGivenNoPred;
            }
        case DIE:
            if (event.agent().type() == PREDATOR) {
                return (event.agent().hasAnySurrounding(PREY, event.time(), trajectory)) ? lpPredDeathGivenPrey : lpPredDeathGivenNoPrey;
            } else {
                return (event.agent().hasAnySurrounding(PREDATOR, event.time(), trajectory)) ? lpPreyDeathGivenPred : lpPreyDeathGivenNoPred;
            }
    }
    assert(false);
    return -INFINITY;
}

template<int GRIDSIZE>
TrajectoryDependencies<PredPreyAgent<GRIDSIZE>> PredPreyAgent<GRIDSIZE>::eventProbDependencies(const Event<PredPreyAgent<GRIDSIZE>> &event) {
    switch(event.act()) {
        case MOVEUP:
        case MOVEDOWN:
        case MOVELEFT:
        case MOVERIGHT:
            return {{},{}};
        case GIVEBIRTH:
        case DIE:
            PredPreyAgent<GRIDSIZE> agent = event.agent();
            Type affectedType = (agent.type()==PREDATOR?PREY:PREDATOR);
            return {{
                            State(event.time(), PredPreyAgent<GRIDSIZE>(agent.xLeft(), agent.yPosition(), affectedType)),
                            State(event.time(), PredPreyAgent<GRIDSIZE>(agent.xRight(), agent.yPosition(), affectedType)),
                            State(event.time(), PredPreyAgent<GRIDSIZE>(agent.xPosition(), agent.yUp(), affectedType)),
                            State(event.time(), PredPreyAgent<GRIDSIZE>(agent.xPosition(), agent.yDown(), affectedType))
                    },
                    {}};
    }
    assert(false);
    return {{},{}};
}


// Should be factored approximation of multinomial given an agent in start position and
// that the constraints are satisfied
// (i.e. given that the prob of this act is not zero)
//template<int GRIDSIZE>
//double PredPreyAgent<GRIDSIZE>::marginalTimestep(Act act) const {
//    constexpr double pPreyNeighbourGivenPred = (1.0-pNoPrey4);
//    constexpr double pPredNeighbourGivenPrey = (1.0-pNoPred4);
//    constexpr double pPreyMove = 0.25 * (1.0 - pPreyDie - pPreyEatenGivenPred*pPredNeighbourGivenPrey - pPreyBirth);
//    constexpr double pPredMove = 0.25 * (1.0 - pPredDie - pPredBirthGivenPrey*pPreyNeighbourGivenPred);
//
//    if(type() == PREDATOR) {
//        switch (act) {
//            case DIE:       return pPredDie;
//            case GIVEBIRTH: return pPredBirthGivenPrey; // assume feasible
//        }
//        return pPredMove;
//    }
//    switch (act) {
//        case DIE:       return pPreyDie + pPreyEatenGivenPred*pPredNeighbourGivenPrey;
//        case GIVEBIRTH: return pPreyBirth;
//    }
//    return pPreyMove;
//}


//template<int GRIDSIZE>
//std::vector<Constraint<ABM::occupation_type>> PredPreyAgent<GRIDSIZE>::constraints(int time, PredPreyAgent<GRIDSIZE>::Act act) const {
//    if(type() == PREDATOR && act == GIVEBIRTH) {
//        return std::vector<Constraint<ABM::occupation_type>>({
//                                                    1*State(time,PredPreyAgent<GRIDSIZE>(xLeft(),yPosition(), PREY)) +
//                                                    1*State(time,PredPreyAgent<GRIDSIZE>(xRight(),yPosition(), PREY)) +
//                                                    1*State(time,PredPreyAgent<GRIDSIZE>(xPosition(),yUp(),PREY)) +
//                                                    1*State(time,PredPreyAgent<GRIDSIZE>(xPosition(),yDown(),PREY)) >= 1
//                                            });
//    }
//    return std::vector<Constraint<ABM::occupation_type>>();
//}

//template<int GRIDSIZE>
//std::vector<PredPreyAgent<GRIDSIZE>> PredPreyAgent<GRIDSIZE>::neighbours() const {
//    Type affectedType = (type()==PREDATOR?PREY:PREDATOR);
//    return {
//        PredPreyAgent<GRIDSIZE>(xLeft(),yPosition(),affectedType),
//        PredPreyAgent<GRIDSIZE>(xRight(),yPosition(),affectedType),
//        PredPreyAgent<GRIDSIZE>(xPosition(),yUp(),affectedType),
//        PredPreyAgent<GRIDSIZE>(xPosition(),yDown(),affectedType)
//    };
//}



// If allowInfeasibleActs is true, then if an act is infeasible then its probability is not set to zero but
// is the expectation averaged over the trajectoryPrior distribution of 'others' (note that this has the
// consequence that the sum of the acts no longer adds to 1 since different acts are relative to different
// distributions of 'others'). This is to enable different uses: for a forward run we set allowInfeasibleActs to
// false and get a LogPMF over acts, whereas when calculating the probability of a exactEndState by setting allowInfeasibleActs
// we extend probability calculation to infeasible trajectories in a reasonably smooth manner.
//template<int GRIDSIZE>
//std::vector<double> PredPreyAgent<GRIDSIZE>::timestep(const ModelState<PredPreyAgent<GRIDSIZE>> &others) const {
//
//    std::vector<double> actDistribution(actDomainSize(),0.0);
//    if (type() == PREDATOR) {
//        if (surroundingCountOf(PREY, others) >= 1) {
//            actDistribution[GIVEBIRTH] = lpPredBirthGivenPrey;
//            actDistribution[DIE] = lpPredDeathGivenPrey;
//        } else {
//            actDistribution[GIVEBIRTH] = lpPredBirthGivenNoPrey;
//            actDistribution[DIE] = lpPredDeathGivenNoPrey;
//        }
//    } else {
//        if (surroundingCountOf(PREDATOR, others) >= 1) {
//            actDistribution[DIE] = lpPreyDeathGivenPred;
//            actDistribution[GIVEBIRTH] = lpPreyBirthGivenPred;
//        } else {
//            actDistribution[DIE] = lpPreyDeathGivenNoPred;
//            actDistribution[GIVEBIRTH] = lpPreyBirthGivenNoPred;
//        }
//    }
//    actDistribution[MOVEUP] = lpMove;
//    actDistribution[MOVEDOWN] = lpMove;
//    actDistribution[MOVELEFT] = lpMove;
//    actDistribution[MOVERIGHT] = lpMove;
//
//
//    assert(fabs(exp(actDistribution[MOVEUP]) + exp(actDistribution[MOVEDOWN]) + exp(actDistribution[MOVELEFT]) + exp(actDistribution[MOVERIGHT]) + exp(actDistribution[DIE]) + exp(actDistribution[GIVEBIRTH]) - 1.0) < 1e-8);
//
////    std::cout << "PPredNeighbour = " << nPredNeighbour*1.0/nPredCalls << " PPreyNeighbour = " << nPreyNeighbour*1.0/nPreyCalls << std::endl;
////    std::cout << "Total agents = " << (others * 1.0).sum() << std::endl;
//    return actDistribution;
//}

#endif //GLPKTEST_PREDPREYAGENT_H
