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

    static constexpr double pNoPred = 0.95; // steady state probability that there are no predators on a square
    static constexpr double pNoPrey = 0.95; // steady state probability that there are no prey on a square
    static constexpr double clustering = 0.3; // tendency to cluster / instability
    static constexpr double pNoPred4 = pNoPred*pNoPred*pNoPred*pNoPred; // no pow function for constexpr :-(
    static constexpr double pNoPrey4 = pNoPrey*pNoPrey*pNoPrey*pNoPrey;

    // values for which there is a steady state prior of uniform, uncorrelated Poisson distribution
    static constexpr double pPredBirthGivenPrey = clustering;
    static constexpr double pPredDie = (1.0-pNoPrey4)*pPredBirthGivenPrey;
    static constexpr double pPreyEatenGivenPred = clustering;
    static constexpr double pPreyDie = 0.1; // stability/mixing?
    static constexpr double pPreyBirth = pPreyDie + (1.0-pNoPred4)*pPreyEatenGivenPred;

    // Agent Domain stuff
    static constexpr int actDomainSize() { return 6; }
};


template<int GRIDSIZE>
class PredPreyAgent: public PredPreyAgentBase {
public:
    // Agent Domain stuff
    static constexpr int domainSize() { return GRIDSIZE*GRIDSIZE*2; }

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
    double marginalTimestep(Act act) const; // returns the probability of this act, given an agent in this state
    std::vector<PredPreyAgent<GRIDSIZE>> consequences(Act act) const; // the consequences of an act
    // returns the constraints implied by the given act of this agent
    std::vector<Constraint<ABM::occupation_type>> constraints(int time, Act act) const; // to be generated automatically by static analysis...eventually.
    std::vector<PredPreyAgent<GRIDSIZE>> neighbours();

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

//    static int nPredCalls = 0;
//    static int nPreyCalls = 0;
//    static int nPredNeighbour = 0;
//    static int nPreyNeighbour = 0;

    std::vector<double> actDistribution(actDomainSize(),0.0);
    if (type() == PREDATOR) {
//        nPredCalls++;
        actDistribution[DIE] = pPredDie;
        if (surroundingCountOf(PREY, others) >= 1) {
//            nPredNeighbour++;
            actDistribution[GIVEBIRTH] = pPredBirthGivenPrey;
        }
    } else {
//        nPreyCalls++;
        actDistribution[GIVEBIRTH] = pPreyBirth;
        actDistribution[DIE] = pPreyDie;
        if (surroundingCountOf(PREDATOR, others) >= 1) {
//            nPreyNeighbour++;
            actDistribution[DIE] += pPreyEatenGivenPred;
        }
    }
    double moveProb = 0.25 * (1.0 - actDistribution[DIE] - actDistribution[GIVEBIRTH]);
    actDistribution[MOVEUP] = moveProb;
    actDistribution[MOVEDOWN] = moveProb;
    actDistribution[MOVELEFT] = moveProb;
    actDistribution[MOVERIGHT] = moveProb;

//    std::cout << "PPredNeighbour = " << nPredNeighbour*1.0/nPredCalls << " PPreyNeighbour = " << nPreyNeighbour*1.0/nPreyCalls << std::endl;
//    std::cout << "Total agents = " << (others * 1.0).sum() << std::endl;
    return actDistribution;
}


// Should be factored approximation of multinomial given an agent in start position and
// that the constraints are satisfied
// (i.e. given that the prob of this act is not zero)
template<int GRIDSIZE>
double PredPreyAgent<GRIDSIZE>::marginalTimestep(Act act) const {
    constexpr double pPreyNeighbourGivenPred = (1.0-pNoPrey4);
    constexpr double pPredNeighbourGivenPrey = (1.0-pNoPred4);
    constexpr double pPreyMove = 0.25 * (1.0 - pPreyDie - pPreyEatenGivenPred*pPredNeighbourGivenPrey - pPreyBirth);
    constexpr double pPredMove = 0.25 * (1.0 - pPredDie - pPredBirthGivenPrey*pPreyNeighbourGivenPred);

    if(type() == PREDATOR) {
        switch (act) {
            case DIE:       return pPredDie;
            case GIVEBIRTH: return pPredBirthGivenPrey; // assume feasible
        }
        return pPredMove;
    }
    switch (act) {
        case DIE:       return pPreyDie + pPreyEatenGivenPred*pPredNeighbourGivenPrey;
        case GIVEBIRTH: return pPreyBirth;
    }
    return pPreyMove;
}


template<int GRIDSIZE>
std::vector<Constraint<ABM::occupation_type>> PredPreyAgent<GRIDSIZE>::constraints(int time, PredPreyAgent<GRIDSIZE>::Act act) const {
    if(type() == PREDATOR && act == GIVEBIRTH) {
        return std::vector<Constraint<ABM::occupation_type>>({
                                                    1*State(time,PredPreyAgent<GRIDSIZE>(xLeft(),yPosition(), PREY)) +
                                                    1*State(time,PredPreyAgent<GRIDSIZE>(xRight(),yPosition(), PREY)) +
                                                    1*State(time,PredPreyAgent<GRIDSIZE>(xPosition(),yUp(),PREY)) +
                                                    1*State(time,PredPreyAgent<GRIDSIZE>(xPosition(),yDown(),PREY)) >= 1
                                            });
    }
    return std::vector<Constraint<ABM::occupation_type>>();
}

template<int GRIDSIZE>
std::vector<PredPreyAgent<GRIDSIZE>> PredPreyAgent<GRIDSIZE>::neighbours() {
    Type affectedType = (type()==PREDATOR?PREY:PREDATOR);
    return {
        PredPreyAgent<GRIDSIZE>(xLeft(),yPosition(),affectedType),
        PredPreyAgent<GRIDSIZE>(xRight(),yPosition(),affectedType),
        PredPreyAgent<GRIDSIZE>(xPosition(),yUp(),affectedType),
        PredPreyAgent<GRIDSIZE>(xPosition(),yDown(),affectedType)
    };
}

#endif //GLPKTEST_PREDPREYAGENT_H
