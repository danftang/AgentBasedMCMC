// Calculates the log of the ratio P(X)/P_i(X) between the real probability and the
// log-linear approximation of the forecast given the start state.
//
// This amounts to \sum_{t \psi a} (\Psi^t_\psi)! log(\pi(\psi, \Psi^t, a)/c^t_{\psi a})
//
// Created by daniel on 15/03/2022.
//

#ifndef ABMCMC_TRAJECTORYIMPORTANCE_H
#define ABMCMC_TRAJECTORYIMPORTANCE_H

#include "Trajectory.h"
#include "IndexSet.h"
#include "ABM.h"
#include "PerturbableFunction.h"

template<class AGENT>
class TrajectoryImportance: public PerturbableFunction<ABM::occupation_type,double> {
public:

    double                       alpha;// = 1.5;// 2.0;  // amplitude of P_i
    StateTrajectory<AGENT>       stateTrajectory;
    std::vector<double>          stateLogProbs;     // factors in the log probability by stateId
    double                       totalLogProb;
    int                          eventTrajecotrySize;

    IndexSet                     staleStates;     // elements of stateTrajectory that need updating, by stateId
    IndexSet                     staleFactors;     // elements of stateLogProbs that need updating, by stateId

    bool                         perturbationsHeldForUndo; //
    std::vector<ABM::occupation_type>   stateOccupationsForUndo;
    std::vector<double>                 factorValuesForUndo;
    double                       totalLogProbForUndo;


//    TrajectoryImportance(const Trajectory<AGENT> &eventTrajectory):
//        eventTrajectory(eventTrajectory),
//        stateTrajectory(eventTrajectory),
//        stateLogProbs(eventTrajectory.nTimesteps()*AGENT::domainSize(),0.0),
//        staleStates(stateLogProbs.size()),
//        staleFactors(stateLogProbs.size()) {
//        initLogProbs();
//    }

    TrajectoryImportance(int nTimesteps, double alpha):
            alpha(alpha),
            stateTrajectory(nTimesteps),
            stateLogProbs(nTimesteps*AGENT::domainSize(),0.0),
            staleStates(stateLogProbs.size()),
            staleFactors(stateLogProbs.size()),
            eventTrajecotrySize(nTimesteps*AGENT::domainSize()*AGENT::actDomainSize()),
            perturbationsHeldForUndo(false) { }


    // register a perturbation to the state
    void perturb(const std::vector<ABM::occupation_type> &X, const std::vector<int> &perturbedIndices) {
//        std::cout << "perturbing" << std::endl;
        if(perturbationsHeldForUndo) clearPerturbations();
        if(hasUndoInformation()) clearUndoInformation();
        for(int i: perturbedIndices) {
            if(i < eventTrajecotrySize) {
                Event<AGENT> changedEvent(i);
                int state = stateId(changedEvent.time(), changedEvent.agent());
                staleStates.insert(state);
                staleFactors.insert(state);
                for (AGENT neighbour: changedEvent.agent().neighbours()) {
                    staleFactors.insert(stateId(changedEvent.time(), neighbour));
                }
            }
        }
    }

    void perturbWithUndo(const std::vector<ABM::occupation_type> &X, const std::vector<int> &perturbedIndices) {
//        std::cout << "perturbing with undo" << std::endl;
        refresh(X);
        if(perturbationsHeldForUndo) {
            clearPerturbations();
            clearUndoInformation();
        }
        totalLogProbForUndo = totalLogProb;
        for(int i : perturbedIndices) {
            if(i < eventTrajecotrySize) {
                Event<AGENT> changedEvent(i);
                int time = changedEvent.time();
                AGENT agent = changedEvent.agent();
                int state = stateId(time, agent);
                if(staleStates.insert(state)) stateOccupationsForUndo.push_back(stateTrajectory[time][agent]);
                if(staleFactors.insert(state)) factorValuesForUndo.push_back(stateLogProbs[state]);
                for (AGENT neighbour: changedEvent.agent().neighbours()) {
                    int neighbourState = stateId(changedEvent.time(), neighbour);
                    if(staleFactors.insert(neighbourState)) factorValuesForUndo.push_back(stateLogProbs[neighbourState]);
                }
            }
        }
    }


    void undoLastPerturbation() {
//        std::cout << "undoing" << std::endl;
        if(staleStates.size() != stateOccupationsForUndo.size() || staleFactors.size() != factorValuesForUndo.size())
            throw("Trying to undo last perturbation, but last perturbation was not undoable");
        // undo state occupation perturbations...
        for(int nnz=0; nnz < staleStates.size(); ++nnz) {
            int stateId = staleStates[nnz];
            int time = stateId / AGENT::domainSize();
            int agent = stateId % AGENT::domainSize();
            stateTrajectory[time][agent] = stateOccupationsForUndo[nnz];
        }
        // undo factor perturbations...
        for(int nnz=0; nnz < staleFactors.size(); ++nnz) {
            stateLogProbs[staleFactors[nnz]] = factorValuesForUndo[nnz];
        }
        totalLogProb = totalLogProbForUndo;
        clearPerturbations();
        clearUndoInformation();
    }

    // log of P(X)/P_i(X)
    double getValue(const std::vector<ABM::occupation_type> &eventTrajectory) {
//        std::cout << "getting value" << std::endl;
        refresh(eventTrajectory);
//        std::cout << "Importance is " << exp(totalLogProb) << std::endl;
//        assert(totalLogProb == 0.0);
        return exp(totalLogProb);
    }


    int nTimesteps() { return stateTrajectory.size(); }


//    const Trajectory<AGENT> &state() { return eventTrajectory; }


    // recalculate stale factors and update total log-prob
    void refresh(const std::vector<ABM::occupation_type> &eventTrajectory) {
//        std::cout << "refreshing" << std::endl;
        if(perturbationsHeldForUndo) return;
//        std::cout << "refreshing, initial importance = " << exp(totalLogProb) << std::endl;
        // first update all state occupation numbers
        for(int stateId: staleStates) {
            int time = stateId / AGENT::domainSize();
            int agent = stateId % AGENT::domainSize();
            stateTrajectory[time][agent] = State<AGENT>(time,agent).forwardOccupation(eventTrajectory);
        }
        // now re-calculate factors
        for(int stateId: staleFactors) {
            int time = stateId / AGENT::domainSize();
            AGENT agent = stateId % AGENT::domainSize();
            int stateOccupation = stateTrajectory[time][agent];
//            std::cout << "state occupation = " << stateOccupation << " time = " << time << " agent = " << agent << std::endl;
            double newLogP = 0.0;
            if(stateOccupation > 0) {
                std::vector<double> actPMF = agent.timestep(stateTrajectory[time]);
                if(stateOccupation > 1) newLogP = lgamma(stateOccupation + 1.0); // log of Phi factorial
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    int actOccupation = eventTrajectory[Event<AGENT>(time, agent, act).id];
                    if(actOccupation > 0 && actPMF[act] > 0.0)
                        newLogP += actOccupation * log(actPMF[act] / agent.marginalTimestep(act));
                }
            }
//            std::cout << "newLogP = " << newLogP << " oldLogP = " << stateLogProbs[stateId] << std::endl;
            assert(newLogP != NAN);
            totalLogProb += newLogP - stateLogProbs[stateId];
            stateLogProbs[stateId] = newLogP;
        }
        if(hasUndoInformation()) {
            perturbationsHeldForUndo = true;
        } else {
            clearPerturbations();
        }
//        std::cout << "Final importance = " << exp(totalLogProb) << std::endl;
    }

    // sets up log probs and state trajectory (non-incremental)
    void setState(const std::vector<ABM::occupation_type> &eventTrajectory) {
        totalLogProb = -log(alpha);
        for (int t = 0; t < nTimesteps(); ++t) {
            for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
                stateTrajectory[t][agentId] = State<AGENT>(t,agentId).forwardOccupation(eventTrajectory);
            }
            for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
                ABM::occupation_type stateOccupation = stateTrajectory[t][agentId];
                double &factorLogProb = stateLogProbs[stateId(t,agentId)];
                factorLogProb = 0.0;
                if(stateOccupation > 0) {
                    AGENT agent(agentId);
                    std::vector<double> actPMF = agent.timestep(stateTrajectory[t]);
                    for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                        double actOccupation = eventTrajectory[Event<AGENT>(t, agentId, actId).id];
                        if (actOccupation > 0 && actPMF[actId] > 0.0) factorLogProb += actOccupation * log(actPMF[actId] / agent.marginalTimestep(actId));
                    }
                }
                if(stateOccupation > 1) factorLogProb += lgamma(stateOccupation + 1.0);
                totalLogProb += factorLogProb;
            }
        }
        stateOccupationsForUndo.clear();
        factorValuesForUndo.clear();
        clearPerturbations();
    }


//    double logPFactor(int factorIndex, const std::vector<double> &state) {
//        int time = factorIndex / AGENT::domainSize();
//        int agent = factorIndex % AGENT::domainSize();
//        std::vector<double> actPMF = AGENT(agent).timestep(stateTrajectory[time]);
//        double logP = lgamma(stateTrajectory[time][agent] + 1.0); // log of Phi factorial
//        for(int act = 0; act < AGENT::actDomainSize(); ++act) {
//            int actOccupation = state[Event<AGENT>(time, agent, act)];
//            if(actOccupation > 0) logP += actOccupation * actPMF[act];
//        }
//        return logP;
//    }


protected:

    void clearPerturbations() {
        staleStates.clear();
        staleFactors.clear();
    }


    void clearUndoInformation() {
        stateOccupationsForUndo.clear();
        factorValuesForUndo.clear();
        perturbationsHeldForUndo = false;
    }

    bool hasUndoInformation() {
        return factorValuesForUndo.size() > 0 || stateOccupationsForUndo.size() > 0;
    }

    int stateId(int time, int agent) { return time * AGENT::domainSize() + agent; }
};


#endif //ABMCMC_TRAJECTORYIMPORTANCE_H
