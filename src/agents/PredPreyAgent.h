//
// Created by daniel on 18/05/2021.
//

// Making AGENTs define their own trajectory accessors?...
//   * need to somehow coordinate with StartStates and Observations
//   * if we define Trajectory[IndexSetFor<Trajectory>] then any linear
//     combination can be extracted from the trajectory and a new index type
//     can be defined as long as it has a cast to IndexSetFor<Trajectory> for
//     each type of trajecotry.
//
//   * alternatively, we go the other way around and define Accessor[trajectory],
//     which probably makes more sense. The accessor can:
//        - calculate occupation number given a trajectory
//        - give it's dependencies for a given trajectory type e.g. dependencies<Trajectory>()
//        - give it's coefficients for a given trajectory type e.g. coefficients<Trajectory>()
//
//   * alternatively, we have accessors in the AGENT definition, so
//        - AGENT::occupation(Trajectory, IndexType)
//        - AGENT::coefficients(Trajectory, IndexType)
//
//  * alternatively, we assume that trajectories have act and state accessors, then let the
//    agent define it's own extra accessors
//
//  * alternatively, we make it easy to add accessors to trajectories...and each AGENT provides
//    its own trajectory type at compile time, and can then assume it's dealing with this
//    trajectory type.
//    A trajectory can be thought odf as trajectory[ACCESSORTYPE]<[time]><[agent]><[act]>
//    or if we assume all accessors are within-timestep then trajectory[time][ACCESSORTYPE]<[agent]><[act]>
//    So we can build a trejectory by supplying a set of accessor types in a pack
//    Each accessor type supplies its own domainSize, a cast to int (indexOf) and linear constraints
//    in terms of other accessors...So accesrots have prerequisites in terms of which
//    Accessors are needed to express their coefficients. So, a State accessor requires
//    an act accessor, a surrounding accessor requires a state accessor
//
//    So... an accessor is a map from various objects to occupation and from int to occupation

#ifndef GLPKTEST_PREDPREYAGENT_H
#define GLPKTEST_PREDPREYAGENT_H

#include <vector>
#include <map>
#include "../ModelState.h"
#include "../version1/Constraint.h"
#include "../ExtendedTrajectory.h"

class PredPreyAgentBase {
public:
    typedef int Act;

    enum ActNames {
        GIVEBIRTH,
        MOVELEFT,
        MOVERIGHT,
        MOVEUP,
        MOVEDOWN,
        DIE
    };
    static constexpr int actDomainSize= 6;


    enum Type {
        PREDATOR,
        PREY
    };
    static constexpr int typeDomainSize = 2;


    static constexpr double pPred = 0.05; // steady state probability that there are no predators on a square
    static constexpr double pPrey = 0.05; // steady state probability that there are no prey on a square
    static constexpr double pNoPred = 1.0-pPred; // steady state probability that there are no predators on a square
    static constexpr double pNoPrey = 1.0-pPrey; // steady state probability that there are no prey on a square
    static constexpr double clustering = 0.33; // tendency to cluster
    static constexpr double pNoPredClose = pNoPred*pNoPred*pNoPred*pNoPred; // no pow function for constexpr :-(
    static constexpr double pNoPreyClose = pNoPrey*pNoPrey*pNoPrey*pNoPrey;
    static constexpr double pPredClose = 1.0 - pNoPredClose;
    static constexpr double pPreyClose = 1.0 - pNoPreyClose;
    static constexpr double APred = 0.666;  // predator birth/death ratio in the presence of prey, 0.5 < APred <= 1.0
    static constexpr double APrey = 0.333; // prey birth/death ratio in the presence of predators, 0.0 <= APrey < 0.5

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

    // marginals
    static const double lpPredBirth;
    static const double lpPreyBirth;
    static const double lpPredDeath;
    static const double lpPreyDeath;

    // Agent Domain stuff
};

// These values ensure that the average population density of pred and prey stays constant in an infinitely large simulation
inline const double PredPreyAgentBase::lpPredBirthGivenPrey = log(clustering*APred);
inline const double PredPreyAgentBase::lpPredDeathGivenPrey = log(clustering*(1.0 - APred));
inline const double PredPreyAgentBase::lpPredBirthGivenNoPrey = log(0.5 * clustering * (1.0 + (1.0-2.0*APred) * pPreyClose / pNoPreyClose));
inline const double PredPreyAgentBase::lpPredDeathGivenNoPrey = log(clustering - exp(lpPredBirthGivenNoPrey));

inline const double PredPreyAgentBase::lpPreyBirthGivenPred = log(clustering*APrey);
inline const double PredPreyAgentBase::lpPreyDeathGivenPred = log(clustering*(1.0-APrey));
inline const double PredPreyAgentBase::lpPreyBirthGivenNoPred = log(0.5 * clustering * (1.0 + (1.0-2.0*APrey) * pPredClose / pNoPredClose));
inline const double PredPreyAgentBase::lpPreyDeathGivenNoPred = log(clustering - exp(lpPreyBirthGivenNoPred));

inline const double PredPreyAgentBase::lpMove = log(0.25 * (1.0 - clustering));

// marginal values
inline const double PredPreyAgentBase::lpPredBirth = log(exp(lpPredBirthGivenPrey)*(1.0-pNoPreyClose) + exp(lpPredBirthGivenNoPrey)*pNoPreyClose);
inline const double PredPreyAgentBase::lpPreyBirth = log(exp(lpPreyBirthGivenPred)*(1.0-pNoPredClose) + exp(lpPreyBirthGivenNoPred)*pNoPredClose);
inline const double PredPreyAgentBase::lpPredDeath = log(exp(lpPredDeathGivenPrey)*(1.0-pNoPreyClose) + exp(lpPredDeathGivenNoPrey)*pNoPreyClose);
inline const double PredPreyAgentBase::lpPreyDeath = log(exp(lpPreyDeathGivenPred)*(1.0-pNoPredClose) + exp(lpPreyDeathGivenNoPred)*pNoPredClose);


template<int GRIDSIZE>
class PredPreyAgent: public PredPreyAgentBase {
public:

    // Agent Domain stuff
    static constexpr int gridsize = GRIDSIZE;
    static constexpr int domainSize = GRIDSIZE*GRIDSIZE*2;

    int stateId;

    PredPreyAgent(int ordinal): stateId(ordinal) {}
    PredPreyAgent( int x, int y, Type type): stateId(x + GRIDSIZE*y + GRIDSIZE*GRIDSIZE*type) { }

//    std::vector<double> timestep(const ModelState<PredPreyAgent<GRIDSIZE>> &others) const;
    std::vector<PredPreyAgent<GRIDSIZE>> consequences(Act act) const; // the consequences of an act

    template<class TRAJECTORY>
    static double logEventProb(const Event<PredPreyAgent<GRIDSIZE>> &event, const TRAJECTORY &trajectory);

//    template<class TRAJECTORY>
//    static std::pair<double,bool> widenedLogEventProb(const Event<PredPreyAgent<GRIDSIZE>> &event, const TRAJECTORY &trajectory);


    template<class DOMAIN>
    static std::vector<int> eventProbDependencies(const Event<PredPreyAgent<GRIDSIZE>> &event);

//    double logMarginalTimestep(Act act) const;

//    template<class TRAJECTORY>
//    bool hasAnySurrounding(Type type, int time, const TRAJECTORY &trajectory) const {
//        return
//                trajectory[State(time,PredPreyAgent(xRight(),yPosition(),type))] ||
//                trajectory[State(time,PredPreyAgent(xLeft(),yPosition(),type))] ||
//                trajectory[State(time,PredPreyAgent(xPosition(),yUp(),type))] ||
//                trajectory[State(time,PredPreyAgent(xPosition(),yDown(),type))];
//    }


    operator int() const { return stateId; }
    int xPosition() const { return stateId%GRIDSIZE; }
    int yPosition() const { return (stateId/GRIDSIZE)%GRIDSIZE; }
    Type type() const { return Type(stateId/(GRIDSIZE*GRIDSIZE)); }
    Type otherType() const { return type()==PREDATOR?PREY:PREDATOR; }

    int xLeft() const { return((stateId+GRIDSIZE-1)%GRIDSIZE); }
    int xRight() const { return((stateId+1)%GRIDSIZE); }
    int yUp() const { return((stateId/GRIDSIZE+1)%GRIDSIZE); }
    int yDown() const { return((stateId/GRIDSIZE + GRIDSIZE-1)%GRIDSIZE); }

    PredPreyAgent<GRIDSIZE> leftOther() const { return PredPreyAgent<GRIDSIZE>(xLeft(),yPosition(),otherType()); }
    PredPreyAgent<GRIDSIZE> rightOther() const { return PredPreyAgent<GRIDSIZE>(xRight(),yPosition(), otherType()); }
    PredPreyAgent<GRIDSIZE> upOther() const { return PredPreyAgent<GRIDSIZE>(xPosition(),yUp(),otherType()); }
    PredPreyAgent<GRIDSIZE> downOther() const { return PredPreyAgent<GRIDSIZE>(xPosition(),yDown(),otherType()); }

    friend std::ostream &operator <<(std::ostream &out, const PredPreyAgent<GRIDSIZE> &agent) {

        out << (agent.type()==PREY?"PREY":"PRED") << ":(" << agent.xPosition() << "," << agent.yPosition() <<")";
        return out;
    }




    // Count of agents of given type north, south, east and west of this agent.
//    double surroundingCountOf(Type type, const ModelState<PredPreyAgent<GRIDSIZE>> &others) const {
//        return others[PredPreyAgent(xRight(),yPosition(),type)] +
//        others[PredPreyAgent(xLeft(),yPosition(),type)] +
//        others[PredPreyAgent(xPosition(),yUp(),type)] +
//        others[PredPreyAgent(xPosition(),yDown(),type)];
//    }



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
    int surroundingCount;
    switch(event.act()) {
        case MOVEUP:
        case MOVEDOWN:
        case MOVELEFT:
        case MOVERIGHT:
            return lpMove;
        case GIVEBIRTH:
            surroundingCount = trajectory.surroundingCountOf(State<PredPreyAgent<GRIDSIZE>>(event.time(), event.agent()));
            if (event.agent().type() == PREDATOR) return surroundingCount >= 1 ? lpPredBirthGivenPrey : lpPredBirthGivenNoPrey;
            return surroundingCount >= 1 ? lpPreyBirthGivenPred : lpPreyBirthGivenNoPred;
        case DIE:
            surroundingCount = trajectory.surroundingCountOf(State<PredPreyAgent<GRIDSIZE>>(event.time(), event.agent()));
            if (event.agent().type() == PREDATOR) return surroundingCount >= 1 ? lpPredDeathGivenPrey : lpPredDeathGivenNoPrey;
            return surroundingCount >= 1 ? lpPreyDeathGivenPred : lpPreyDeathGivenNoPred;
    }
    assert(false);
    return -INFINITY;
}


//template<int GRIDSIZE>
//template<class TRAJECTORY>
//std::pair<double,bool> PredPreyAgent<GRIDSIZE>::widenedLogEventProb(const Event<PredPreyAgent<GRIDSIZE>> &event, const TRAJECTORY &trajectory) {
//    int surroundingCount;
//    switch(event.act()) {
//        case MOVEUP:
//        case MOVEDOWN:
//        case MOVELEFT:
//        case MOVERIGHT:
//            return std::pair(lpMove,true);
//        case GIVEBIRTH:
//            surroundingCount = trajectory.surroundingCountOf(State<PredPreyAgent<GRIDSIZE>>(event.time(), event.agent()));
//            if (event.agent().type() == PREDATOR) {
////                return (surroundingCount > 0) ? std::pair(lpPredBirthGivenPrey,true):std::pair(lpPredBirthGivenNoPrey,true);
//                return std::pair(lpPredBirth,true);
//            }
////            return (surroundingCount > 0) ? std::pair(lpPreyBirthGivenPred,true) : std::pair(lpPreyBirthGivenNoPred,true);
//            return std::pair(lpPreyBirth,true);
//        case DIE:
//            surroundingCount = trajectory.surroundingCountOf(State<PredPreyAgent<GRIDSIZE>>(event.time(), event.agent()));
//            if (event.agent().type() == PREDATOR) {
////                return (surroundingCount > 0) ? std::pair(lpPredDeathGivenPrey,true) : std::pair(lpPredDeathGivenNoPrey,true);
//                return std::pair(lpPredDeath,true);
//            }
////            return (surroundingCount > 0) ? std::pair(lpPreyDeathGivenPred,true) : std::pair(lpPreyDeathGivenNoPred,true);
//            return std::pair(lpPreyDeath,true);
//    }
//    assert(false);
//    return std::pair(-INFINITY,false);
//}


template<int GRIDSIZE>
template<class DOMAIN>
std::vector<int> PredPreyAgent<GRIDSIZE>::eventProbDependencies(const Event<PredPreyAgent<GRIDSIZE>> &event) {
    switch(event.act()) {
        case MOVEUP:
        case MOVEDOWN:
        case MOVELEFT:
        case MOVERIGHT:
            return {};
        case GIVEBIRTH:
        case DIE:
            return { DOMAIN::surroundingCountIndexOf(State<PredPreyAgent<GRIDSIZE>>(event.time(), event.agent())) };
//            return {}; // TODO: test!!!!
    }
    assert(false);
    return {};
}


// Should be factored approximation of multinomial given an agent in start position and
// that the constraints are satisfied
// (i.e. given that the prob of this act is not zero)
//template<int GRIDSIZE>
//double PredPreyAgent<GRIDSIZE>::logMarginalTimestep(Act act) const {
//    switch(act) {
//        case MOVEUP:
//        case MOVEDOWN:
//        case MOVELEFT:
//        case MOVERIGHT:
//            return lpMove;
//        case GIVEBIRTH:
//            return type() == PREDATOR ? lpPredBirthGivenPrey:lpPreyBirthGivenNoPred;
//        case DIE:
//            return type() == PREDATOR ? lpPredDeathGivenNoPrey:lpPreyDeathGivenPred;
//    }
//    assert(false);
//    return -INFINITY;
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
//    std::vector<double> actDistribution(actDomainSize,0.0);
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
