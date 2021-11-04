//
// Created by daniel on 18/05/2021.
//

#ifndef GLPKTEST_PREDPREYAGENT_H
#define GLPKTEST_PREDPREYAGENT_H

#include <vector>
#include <map>
#include "glpkpp.h"
#include "../ModelState.h"

template<int GRIDSIZE>
class PredPreyAgent {
public:

    typedef int Act;


    enum ActNames {
        DIE,
        MOVELEFT,
        MOVERIGHT,
        MOVEUP,
        MOVEDOWN,
        GIVEBIRTH
    };


    enum Type {
        PREDATOR,
        PREY
    };

    static constexpr double pPredBirthGivenPrey = 0.1; // 0.5; // birth prob given prey
    static constexpr double pPredDie = 0.07; // death prob
    static constexpr double pPreyBirth = 0.06; // birth prob
    static constexpr double pPreyDie = 0.02;//0.03; // death prob
    static constexpr double pPreyEatenGivenPred = 0.15; // 0.55; // death prob given pred


    // Agent Domain stuff
    static int domainSize() { return GRIDSIZE*GRIDSIZE*2; }
    static constexpr int actDomainSize() { return 6; }

    int stateId;

    PredPreyAgent(int ordinal): stateId(ordinal) {}
    PredPreyAgent( int x, int y, Type type): stateId(x + GRIDSIZE*y + GRIDSIZE*GRIDSIZE*type) { }
    operator int() const { return stateId; }
    int xPosition() const { return stateId%GRIDSIZE; }
    int yPosition() const { return (stateId/GRIDSIZE)%GRIDSIZE; }
    Type type() const { return Type(stateId/(GRIDSIZE*GRIDSIZE)); }

    int xLeft() const { return((stateId+GRIDSIZE-1)%GRIDSIZE); }
    int xRight() const { return((stateId+1)%GRIDSIZE); }
    int yUp() const { return((stateId/GRIDSIZE+1)%GRIDSIZE); }
    int yDown() const { return((stateId/GRIDSIZE + GRIDSIZE-1)%GRIDSIZE); }


    std::vector<double> timestep(const ModelState<PredPreyAgent<GRIDSIZE>> &others) const;
    std::vector<double> marginalTimestep() const;
    std::vector<PredPreyAgent<GRIDSIZE>> consequences(Act act) const; // the consequences of an act
    // returns the constraints implied by the given act
    std::vector<glp::Constraint> constraints(int time, Act act) const; // to be generated automatically by static analysis...eventually.

    friend std::ostream &operator <<(std::ostream &out, const PredPreyAgent<GRIDSIZE> &agent) {

        out << (agent.type()==PREY?"PREY":"PRED") << ":(" << agent.xPosition() << "," << agent.yPosition() <<")";
        return out;
    }


    double surroundingCountOf(Type type, const ModelState<PredPreyAgent<GRIDSIZE>> &others) const {
        return others[PredPreyAgent(xRight(),yPosition(),type)] +
        others[PredPreyAgent(xLeft(),yPosition(),type)] +
        others[PredPreyAgent(xPosition(),yUp(),type)] +
        others[PredPreyAgent(xPosition(),yDown(),type)];
    }

//
//
//    // returns a constraint in terms of state occupation numbers
//    // Manually generated for now...


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


// If allowInfeasibleActs is true, then if an act is infeasible then its probability is not set to zero but
// is the expectation averaged over the trajectoryPrior distribution of 'others' (note that this has the
// consequence that the sum of the acts no longer adds to 1 since different acts are relative to different
// distributions of 'others'). This is to enable different uses: for a forward run we set allowInfeasibleActs to
// false and get a LogPMF over acts, whereas when calculating the probability of a exactEndState by setting allowInfeasibleActs
// we extend probability calculation to infeasible trajectories in a reasonably smooth manner.
template<int GRIDSIZE>
std::vector<double> PredPreyAgent<GRIDSIZE>::timestep(const ModelState<PredPreyAgent<GRIDSIZE>> &others) const {

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

    return actDistribution;
}


template<int GRIDSIZE>
std::vector<double> PredPreyAgent<GRIDSIZE>::marginalTimestep() const {
    constexpr double pPreyNeighbourGivenPred = 0.2;
    constexpr double pPredNeighbourGivenPrey = 0.1;

    std::vector<double> actDistribution(actDomainSize(),0.0);
    if (type() == PREDATOR) {
        actDistribution[DIE] = pPredDie;
        actDistribution[GIVEBIRTH] = pPredBirthGivenPrey*pPreyNeighbourGivenPred;
    } else {
        actDistribution[GIVEBIRTH] = pPreyBirth;
        actDistribution[DIE] = pPreyDie + pPreyEatenGivenPred*pPredNeighbourGivenPrey;
    }
    double moveProb = 0.25 * (1.0 - actDistribution[DIE] - actDistribution[GIVEBIRTH]);
    actDistribution[MOVEUP] = moveProb;
    actDistribution[MOVEDOWN] = moveProb;
    actDistribution[MOVELEFT] = moveProb;
    actDistribution[MOVERIGHT] = moveProb;

    return actDistribution;
}


template<int GRIDSIZE>
std::vector<glp::Constraint> PredPreyAgent<GRIDSIZE>::constraints(int time, PredPreyAgent<GRIDSIZE>::Act act) const {
    if(type() == PREDATOR && act == GIVEBIRTH) {
        return std::vector<glp::Constraint>({
                                                    1.0*State(time,PredPreyAgent<GRIDSIZE>(xLeft(),yPosition(), PREY)) +
                                                    1.0*State(time,PredPreyAgent<GRIDSIZE>(xRight(),yPosition(), PREY)) +
                                                    1.0*State(time,PredPreyAgent<GRIDSIZE>(xPosition(),yUp(),PREY)) +
                                                    1.0*State(time,PredPreyAgent<GRIDSIZE>(xPosition(),yDown(),PREY)) >= 1.0
                                            });
    }
    return std::vector<glp::Constraint>();
}

#endif //GLPKTEST_PREDPREYAGENT_H
