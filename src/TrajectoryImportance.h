//
// Created by daniel on 09/04/2022.
//

#ifndef ABMCMC_TRAJECTORYIMPORTANCE_H
#define ABMCMC_TRAJECTORYIMPORTANCE_H

#include <vector>
#include "StateTrajectory.h"
#include "IndexSet.h"
#include "PerturbableFunction.h"

template<class AGENT>
class TrajectoryImportance: public PerturbableFunction<ABM::occupation_type, double> {
public:

    StateTrajectory<AGENT>       stateTrajectory;       // current state values
    std::vector<double>          stateLogProbs;         // current value of factors in the log probability by stateId()
    double                       totalLogProb;          // current log prob

    IndexSet                     staleStates;     // elements of stateTrajectory that need updating or undoing, by stateId
    IndexSet                     staleFactors;     // elements of stateLogProbs that need updating or undoing, by stateId

    std::vector<ABM::occupation_type>   stateOccupationsForUndo;
    std::vector<double>                 factorValuesForUndo;
    double                              totalLogProbForUndo;


    TrajectoryImportance(int nTimesteps):
            stateTrajectory(nTimesteps),
            stateLogProbs(nTimesteps*AGENT::domainSize(),0.0),
            staleStates(stateLogProbs.size()),
            staleFactors(stateLogProbs.size()) {}


    void perturb(const std::vector<ABM::occupation_type> &X, const std::vector<int> &perturbedIndices) {
//        std::cout << "perturbing with undo" << std::endl;
        clearUndoInformation();

        totalLogProbForUndo = totalLogProb;
        int eventTrajectorySize = nTimesteps()*AGENT::domainSize()*AGENT::actDomainSize();
        for(int i : perturbedIndices) {
            if(i < eventTrajectorySize) {
                Event<AGENT> changedEvent(i);
                State<AGENT> state(changedEvent.time(), changedEvent.agent());
                int stateId = state.id();
                if(staleStates.insert(stateId)) stateOccupationsForUndo.push_back(stateTrajectory[state.time][state.agent]);
                if(staleFactors.insert(stateId)) factorValuesForUndo.push_back(stateLogProbs[stateId]);
                for (AGENT neighbour: changedEvent.agent().neighbours()) {
                    int neighbourState = State<AGENT>(state.time, neighbour).id();
                    if(staleFactors.insert(neighbourState)) factorValuesForUndo.push_back(stateLogProbs[neighbourState]);
                }
            }
        }

        // now populate caches
        stateOccupationsForUndo.reserve(staleStates.size());
        for(int stateId: staleStates) {
            State<AGENT> state(stateId);
            stateOccupationsForUndo.push_back(stateTrajectory[state.time][state.agent]);
        }
        factorValuesForUndo.reserve(staleFactors.size());
        for(int stateId: staleFactors) {
            stateOccupationsForUndo.push_back(stateLogProbs[stateId]);
        }

    }


    void undo() {
//        std::cout << "undoing" << std::endl;
        for(int nnz=0; nnz < staleStates.size(); ++nnz) {
            State<AGENT> state(staleStates[nnz]);
            stateTrajectory[state.time][state.agent] = stateOccupationsForUndo[nnz];
        }
        for(int nnz=0; nnz < staleFactors.size(); ++nnz) {
            stateLogProbs[staleFactors[nnz]] = factorValuesForUndo[nnz];
        }
        totalLogProb = totalLogProbForUndo;
        clearUndoInformation();
    }

    // log of P(X)/P_i(X)
    double getLogValue(const std::vector<ABM::occupation_type> &eventTrajectory) {
//        std::cout << "getting value" << std::endl;
        refresh(eventTrajectory);
//        std::cout << "Importance is " << exp(totalLogProb) << std::endl;
//        assert(totalLogProb == 0.0);
//        sanityCheck(eventTrajectory);
        return totalLogProb;
    }


    int nTimesteps() { return stateTrajectory.size(); }


//    const Trajectory<AGENT> &state() { return eventTrajectory; }


    // recalculate stale factors and update total log-prob
    void refresh(const std::vector<ABM::occupation_type> &eventTrajectory) {
//        std::cout << "refreshing, initial importance = " << exp(totalLogProb) << std::endl;
        // first update all state occupation numbers
        for(int stateId: staleStates) {
            State<AGENT> state(stateId);
//            std::cout << "Refreshing state occupation " << State<AGENT>(time,agent) << std::endl;
            stateTrajectory[state.time][state.agent] = state.fermionicBoundedForwardOccupation(eventTrajectory);
        }
        // now re-calculate factors
        for(int stateId: staleFactors) {
            State<AGENT> state(stateId);
//            std::cout << "Refreshing factor " << State<AGENT>(time,agent) << std::endl;
            int stateOccupation = stateTrajectory[state.time][state.agent];
//            std::cout << "state occupation = " << stateOccupation << " time = " << time << " agent = " << agent << std::endl;
            double newLogP = (stateOccupation > 1)?lgamma(stateOccupation + 1.0):0.0; // log of Phi factorial
            if(stateOccupation > 0) {
                std::vector<double> actPMF = state.agent.timestep(stateTrajectory[state.time]);
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    int actOccupation = eventTrajectory[Event<AGENT>(state.time, state.agent, act).id];
                    if(actOccupation > 0 && actPMF[act] > 0.0) {
//                        std::cout << "agent =" << agent << " act =" << act << " P=" << actPMF[act] << " P_i=" << agent.marginalTimestep(act) << " P/P_i = " << actPMF[act] / agent.marginalTimestep(act) << std::endl;
//                        newLogP += actOccupation * log(actPMF[act] / agent.marginalTimestep(act));
                        newLogP += log(actPMF[act] / state.agent.marginalTimestep(act)); // Fermionic bounding
                    }
                }
            }
//            std::cout << "newLogP = " << newLogP << " oldLogP = " << stateLogProbs[stateId] << std::endl;
            totalLogProb += newLogP - stateLogProbs[stateId];
            stateLogProbs[stateId] = newLogP;
        }
//        std::cout << "Final importance = " << exp(totalLogProb) << std::endl;
    }

    // sets up log probs and state trajectory (non-incremental)
    void setState(const std::vector<ABM::occupation_type> &eventTrajectory) {
        totalLogProb = 0.0;//-log(alpha);
        for (int t = 0; t < nTimesteps(); ++t) {
            for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
                stateTrajectory[t][agentId] = State<AGENT>(t,agentId).fermionicBoundedForwardOccupation(eventTrajectory);
            }
            for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
                ABM::occupation_type stateOccupation = stateTrajectory[t][agentId];
                if(stateOccupation > 0) {
                    State<AGENT> state(t,agentId);
                    double &factorLogProb = stateLogProbs[state.id()];
                    factorLogProb = (stateOccupation > 1)?lgamma(stateOccupation + 1.0):0.0;
                    std::vector<double> actPMF = state.agent.timestep(stateTrajectory[t]);
                    for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                        double actOccupation = eventTrajectory[Event<AGENT>(t, agentId, actId).id];
                        if (actOccupation > 0 && actPMF[actId] > 0.0) {
                            factorLogProb += log(actPMF[actId] / state.agent.marginalTimestep(actId)); // Fermionic bounding
                        }
                    }
                    totalLogProb += factorLogProb;
                }
            }
        }
        clearUndoInformation();
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


    void clearUndoInformation() {
        staleStates.clear();
        staleFactors.clear();
        stateOccupationsForUndo.clear();
        factorValuesForUndo.clear();
    }

//    bool hasUndoInformation() {
//        return factorValuesForUndo.size() > 0 || stateOccupationsForUndo.size() > 0;
//    }

//    int stateId(int time, int agent) { return time * AGENT::domainSize() + agent; }

    void sanityCheck(const std::vector<ABM::occupation_type> &eventTrajectory) {
        double logP = 0.0;//-log(alpha);
        int tEnd = nTimesteps();
        for (int t = 0; t < tEnd; ++t) {
            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId)
                assert(stateTrajectory[t][agentId] == State<AGENT>(t,agentId).fermionicBoundedForwardOccupation(eventTrajectory));
            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                ABM::occupation_type stateOccupation = stateTrajectory[t][agentId];
                if(stateOccupation > 0) {
                    AGENT agent(agentId);
                    std::vector<double> actPMF = agent.timestep(stateTrajectory[t]);
                    double agentStateLogP = (stateOccupation > 1)?lgamma(stateOccupation + 1.0):0.0;
                    for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                        double actOccupation = eventTrajectory[Event<AGENT>(t, agentId, actId).id];
                        if (actOccupation > 0 && actPMF[actId] > 0.0)
                            agentStateLogP += log(actPMF[actId] / agent.marginalTimestep(actId)); // Fermionic bounding
                    }
                    assert(fabs(agentStateLogP - stateLogProbs[State<AGENT>(t,agentId).id()]) < 1e-5);
                    logP += agentStateLogP;
                }
            }
        }

        if(fabs(logP - totalLogProb) > 1e-5) {
            std::cout << "Failed totalLogProb sanity check. True logP = " << logP << " cached LogP = " << totalLogProb <<std::endl;
        }
//        std::cout << "Sanity check logP=" << logP << " cahchedLogProb=" << totalLogProb << std::endl;
        assert(fabs(logP-totalLogProb) < 1e-5);
    }
};


#endif //ABMCMC_TRAJECTORYIMPORTANCE_H
