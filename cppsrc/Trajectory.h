//
// Created by daniel on 04/05/2021.
//

#ifndef GLPKTEST_TRAJECTORY_H
#define GLPKTEST_TRAJECTORY_H

#include "glpkpp.h"
#include "ModelState.h"
#include "StateTrajectory.h"
#include "Random.h"
#include "ActFermionicDistribution.h"

template<typename AGENT>
class Trajectory: public std::vector<double> {
public:
//    static constexpr double tol = SimplexMCMC::tol;

    explicit Trajectory(int nTimesteps): std::vector<double>(nTimesteps*AGENT::domainSize()*AGENT::actDomainSize()+1) { }
    Trajectory(std::vector<double> &&rvalue): std::vector<double>(rvalue) {
        assert((size()-1)%(AGENT::domainSize()*AGENT::actDomainSize()) == 0);
    }
    explicit Trajectory(const std::vector<double> &lvalue): std::vector<double>(lvalue) {
        assert((size()-1)%(AGENT::domainSize()*AGENT::actDomainSize()) == 0);
    }

    // execute forward from a given start state distribution, choosing a solution with probability
    // proportional to the joint times the probability of the start state
    Trajectory(int nTimesteps, const ModelState<AGENT> &startState) : Trajectory(nTimesteps) {
        bool isValid;
        int nAttempts=0;
        do {
            ModelState<AGENT> t0State = startState;
            ModelState<AGENT> t1State;
            isValid = true;
            for (int t = 0; t < nTimesteps; ++t) {
                for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                    AGENT agent(agentId);
                    int nAgents = t0State[agentId];
                    ActFermionicDistribution actPMF = ActFermionicDistribution(agent.timestep(t0State, 0.0));
                    std::vector<bool> acts = actPMF.sampleUnordered(nAgents);
                    if(acts.size() == AGENT::actDomainSize()) {
                        for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                            if (acts[actId]) {
                                operator[](Event(t, agent, actId)) = 1.0;
                                t1State += agent.consequences(actId);
                            } else {
                                operator[](Event(t, agent, actId)) = 0.0;
                            }
                        }
                    } else {
                        isValid = false;
                        t = nTimesteps;
                        agentId = AGENT::domainSize();
                    }
                }
                t0State.setToZero();
                t0State.swap(t1State);
            }
            if(++nAttempts > 1000) throw(std::runtime_error("Can't create act-Fermionic prior sample of Trajectory. Probably too many agents for Fermionicity."));
        } while(!isValid);
    }


    Trajectory(int nTimesteps, const std::function<std::vector<double>()> &startStateSampler) : Trajectory(nTimesteps) {
        bool isValid;
        int nAttempts = 0;
        do {
            ModelState<AGENT> t0State = startStateSampler();
            ModelState<AGENT> t1State;
            isValid = true;
            for (int t = 0; t < nTimesteps; ++t) {
                for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                    AGENT agent(agentId);
                    int nAgents = t0State[agentId];
                    for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                        (*this)[Event(t, agent, actId)] = 0.0;
                    }
                    // now choose acts for each of nAgents
                    std::vector<double> actPMF = agent.timestep(t0State, 0.0);
                    std::discrete_distribution<int> actDistribution(actPMF.begin(), actPMF.end());
                    for(int m=0; m <nAgents; ++m) {
                        int chosenAct = actDistribution(Random::gen);
                        int event = Event(t, agent, chosenAct);
                        if((*this)[event] == 0.0) {
                            (*this)[event] = 1.0;
                            t1State += agent.consequences(chosenAct);
                        } else { // reject
                            isValid = false;
                            m = nAgents;
                            agentId = AGENT::domainSize();
                            t = nTimesteps;
                        }
                    }
                }
                t0State.setToZero();
                t0State.swap(t1State);
            }
            if (++nAttempts > 1000)
                throw (std::runtime_error(
                        "Can't create act-Fermionic prior sample of Trajectory. Too many agents for Fermionicity to be a good assumption."));
        } while (!isValid);
    }

    // event count (use Event(time,agent,act) instead)
//    double operator()(int time, const AGENT &agent, const typename AGENT::Act &act) const {
//        return (*this)[Event(time,agent,act)];
//    }

    // occupation number of an agent state at a particular time
    // time must be between 0 and nTimesteps-1 (final state not implemented at present)
    // use pre-multiply by agent state? (vector dot product)
    double operator()(int time, const AGENT &agent) const {
        int beginIndex = Event(time,agent,0);
        assert(beginIndex < size()); // TODO: Implement final state occupation numbers
        int endIndex = beginIndex + AGENT::actDomainSize();
        double occupation = 0.0;
        int eventId;
        for(int eventId=beginIndex; eventId<endIndex; ++eventId) {
            occupation += (*this)[eventId];
        }
        return occupation;
    }

    // time slice for time is in [0...nTimesteps]
    // i.e. final state is implemented
    // use ModelState constructor instead?
    ModelState<AGENT> operator()(int time) const {
        assert(time>=0 && time<=nTimesteps());
        if(time == nTimesteps()) return endState();
        ModelState<AGENT> timeslice;
        int beginIndex = Event(time,AGENT(0),0);
        int endIndex = beginIndex + AGENT::actDomainSize() * AGENT::domainSize();
        for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
            if (double occupation = (*this)[eventId]; fabs(occupation) > tol)
                timeslice[Event<AGENT>(eventId).agent()] += occupation;
        }
        return timeslice;
    }


    ModelState<AGENT> endState() const {
        ModelState<AGENT> timeslice;
        int beginIndex = Event(nTimesteps()-1,AGENT(0),0);
        int endIndex = beginIndex + AGENT::actDomainSize() * AGENT::domainSize();
        for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
            if (double occupation = (*this)[eventId]; fabs(occupation) > tol) {
                for(const AGENT &consequence : Event<AGENT>(eventId).consequences()) {
                    timeslice[consequence] += occupation;
                }
            }
        }
        return timeslice;
    }


    int nTimesteps() const { return (size()-1)/(AGENT::domainSize()*AGENT::actDomainSize()); }

    double logProb() const {
        const double infeasibilityPenalty = 0.01;

        StateTrajectory<AGENT> stateTrajectory(*this);
        double logP = 0.0;
        int tEnd = nTimesteps();
        for (int t = 0; t < tEnd; ++t) {
            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                std::vector<double> actPMF = AGENT(agentId).timestep(stateTrajectory[t], infeasibilityPenalty);
                for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                    double occupation = fabs((*this)[Event<AGENT>(t, agentId, actId)]);
                    if (occupation > tol) {
                        logP += occupation * log(actPMF[actId]);
//                   - CombinatoricsUtils.factorialLog(occupation) // add this for non-act-fermionic trajectories
                    }
                }
            }
        }

        // If any state occupation number, m, is greater than 1 then we need to
        // multiply the prob by !m since the m agents can be assigned to m acts in
        // !m ways.
        for(const ModelState<AGENT> &step: stateTrajectory) {
            for(double occupation: step) {
                if(fabs(occupation) > 1.0 + tol) {
                    logP += lgamma(fabs(occupation) + 1.0);
                }
            }
        }
//        std::cout << "Trajectory prior = " << logP << std::endl;
        return logP;
    }



//    friend std::ostream &operator <<(std::ostream &out, const Trajectory<AGENT> &solution) {
//        Trajectory sortedTrajectory(solution);
//        sortedTrajectory.sort();
//
//        int timestep = -1;
//        for(auto eventCount: sortedTrajectory) {
//            Event<AGENT> event = eventCount.index;
//            if(event.time() != timestep) {
//                out << "time = " << event.time() << std::endl;
//                timestep = event.time();
//            }
//            out << "  " << eventCount.value << "*" << event << std::endl;
//        }
//        return out;
//    }


};

#endif //GLPKTEST_TRAJECTORY_H
