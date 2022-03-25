//
// Created by daniel on 04/05/2021.
//

#ifndef GLPKTEST_TRAJECTORY_H
#define GLPKTEST_TRAJECTORY_H

#include <cfloat>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include "ModelState.h"
#include "StateTrajectory.h"
#include "Random.h"
#include "ActFermionicDistribution.h"
#include "Sampler.h"
#include "ABM.h"

template<typename AGENT>
class Trajectory: public std::vector<ABM::occupation_type> {
public:
//    static constexpr double tol = SimplexMCMC::tol;

    explicit Trajectory(int nTimesteps): std::vector<value_type>(dimension(nTimesteps)) { }

    explicit Trajectory(std::vector<value_type> &&rvalue): std::vector<value_type>(rvalue) {
        assert((this->size()-1)%(AGENT::domainSize()*AGENT::actDomainSize()) == 0);
    }

    explicit Trajectory(const std::vector<value_type> &lvalue): std::vector<value_type>(lvalue) {
        assert((this->size()-1)%(AGENT::domainSize()*AGENT::actDomainSize()) == 0);
    }

    // execute forward from a given start state distribution, choosing a exactEndState with probability
    // proportional to the joint times the probability of the start state
//    Trajectory(int nTimesteps, const ModelState<AGENT> &startState) : Trajectory(nTimesteps) {
//        bool isValid;
//        int nAttempts=0;
//        do {
//            ModelState<AGENT> t0State = startState;
//            ModelState<AGENT> t1State;
//            isValid = true;
//            for (int t = 0; t < nTimesteps; ++t) {
//                for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
//                    AGENT agent(agentId);
//                    int nAgents = t0State[agentId];
//                    ActFermionicDistribution actPMF = ActFermionicDistribution(agent.timestep(t0State, 0.0));
//                    std::vector<bool> acts = actPMF.sampleUnordered(nAgents);
//                    if(acts.size() == AGENT::actDomainSize()) {
//                        for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
//                            if (acts[actId]) {
//                                operator[](Event(t, agent, actId)) = 1.0;
//                                t1State += agent.consequences(actId);
//                            } else {
//                                operator[](Event(t, agent, actId)) = 0.0;
//                            }
//                        }
//                    } else {
//                        isValid = false;
//                        t = nTimesteps;
//                        agentId = AGENT::domainSize();
//                    }
//                }
//                t0State.setToZero();
//                t0State.swap(t1State);
//            }
//            if(++nAttempts > 1000) throw(std::runtime_error("Can't create act-Fermionic trajectoryPrior sample of Trajectory. Probably too many agents for Fermionicity."));
//        } while(!isValid);
//    }


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
                    // now choose acts for each of nAgents from act Fermionic distribution

                    std::vector<double> actPMF = agent.timestep(t0State);
                    ActFermionicDistribution actDistribution(actPMF);
                    assert(nAgents <= actPMF.size());
                    std::vector<bool> chosenActs = actDistribution.sampleUnordered(nAgents);
                    for(int act=0; act < chosenActs.size(); ++act) {
                        if(chosenActs[act]) {
                            int event = Event(t, agent, act);
                            (*this)[event] = 1.0;
                            t1State += agent.consequences(act);
                        }
                    }

//                    std::vector<double> actPMF = agent.timestep(t0State);
//                    std::discrete_distribution<int> actDistribution(actPMF.begin(), actPMF.end());
//                    for(int m=0; m <nAgents; ++m) { // TODO: replace with ActFermionicDistribution
//                        int chosenAct = actDistribution(Random::gen);
//                        int event = Event(t, agent, chosenAct);
//                        if((*this)[event] == 0.0) {
//                            (*this)[event] = 1.0;
//                            assert(event < size());
//                            t1State += agent.consequences(chosenAct);
//                        } else { // reject
//                            isValid = false;
//                            m = nAgents;
//                            agentId = AGENT::domainSize();
//                            t = nTimesteps;
//                        }
//                    }

                }
                t0State.setToZero();
                t0State.swap(t1State);
            }
            if (++nAttempts > 4000)
                throw (std::runtime_error(
                        "Can't create act-Fermionic trajectoryPrior sample of Trajectory. Too many agents for Fermionicity to be a good assumption."));
        } while (!isValid);
    }

    // event count (use Event(time,agent,act) instead)
//    double operator()(int time, const AGENT &agent, const typename AGENT::Act &act) const {
//        return (*this)[Event(time,agent,act)];
//    }

    // occupation number of an agent state at a particular time
    value_type operator[](const State<AGENT> &state) const {
        assert(state.time >= 0 && state.time <= nTimesteps());
        return(time < nTimesteps())?state.forwardOccupation(*this):state.backwardOccupation(*this);
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
            if (value_type occupation = (*this)[eventId]; occupation > 0)
                timeslice[Event<AGENT>(eventId).agent()] += occupation;
        }
        return timeslice;
    }


    ModelState<AGENT> endState() const {
        ModelState<AGENT> timeslice;
        int beginIndex = Event(nTimesteps()-1,AGENT(0),0);
        int endIndex = beginIndex + AGENT::actDomainSize() * AGENT::domainSize();
        for (int eventId = beginIndex; eventId < endIndex; ++eventId) {
            if (value_type occupation = (*this)[eventId]; occupation != 0) {
                for(const AGENT &consequence : Event<AGENT>(eventId).consequences()) {
                    timeslice[consequence] += occupation;
                }
            }
        }
        return timeslice;
    }


    value_type forwardOccupationNumber(const State<AGENT> &agentState) const;

    value_type backwardOccupationNumber(const State<AGENT> &agentState) const;

    value_type occupationNumber(const State<AGENT> &agentState) const;


    int nTimesteps() const { return (size()-1)/(AGENT::domainSize()*AGENT::actDomainSize()); }

    int dimension() const { return size(); }

    static int dimension(int nTimesteps) { return AGENT::domainSize()*AGENT::actDomainSize()*nTimesteps + 1; }

    // Log probability given fixed start state, as defined by this trajectory
//    double logProbOld() const {
////        const double infeasibilityPenalty = 0.5;
//
//        StateTrajectory<AGENT> stateTrajectory(*this);
//        double logP = 0.0;
//        int tEnd = nTimesteps();
//        for (int t = 0; t < tEnd; ++t) {
//            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
//                std::vector<double> actPMF = AGENT(agentId).timestep(stateTrajectory[t]);
//                double agentStateLogP;
//                do {
//                    agentStateLogP = 0.0;
//                    for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
////                        double occupation = fabs((*this)[Event<AGENT>(t, agentId, actId)]);
//                        double occupation = (*this)[Event<AGENT>(t, agentId, actId)];
//                        if(occupation < 0.0) occupation = 0.0; else if(occupation > 1.0) occupation = 1.0; // Clamp to Fermionic limits
//                        if (occupation > tol) {
//                            agentStateLogP += occupation * log(actPMF[actId]);
//                            // - CombinatoricsUtils.factorialLog(occupation) // add this for non-act-fermionic trajectories
//                        }
//                    }
//                    if(agentStateLogP <= -DBL_MAX) actPMF = AGENT(agentId).marginalTimestep();
//                } while(agentStateLogP <= -DBL_MAX);
//                logP += agentStateLogP;
//            }
//        }
//
//        // If any state occupation number, m, is greater than 1 then we need to
//        // multiply the prob by !m since the m agents can be assigned to m acts in
//        // !m ways.
//        for(const ModelState<AGENT> &step: stateTrajectory) {
//            for(double occupation: step) {
//                if(fabs(occupation) > 1.0 + tol) {
//                    logP += lgamma(fabs(occupation) + 1.0);
//                }
//            }
//        }
////        std::cout << "Trajectory trajectoryPrior = " << logP << std::endl;
//        return logP;
//    }


    // Generate random set of observations of this trajectory
    std::vector<AgentStateObservation<AGENT>> generateObservations(double pMakeObservation, double pObserveIfPresent) {
        std::vector<AgentStateObservation<AGENT>> observations;
        int T = nTimesteps();
        for (int t=0; t<T;++t) {
            for (int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
                if (Random::nextDouble() < pMakeObservation) {
                    AGENT agent(agentId);
                    int nObserved = Random::nextBinomial((*this)[State<AGENT>(t,agent)], pObserveIfPresent);
                    observations.emplace_back(State<AGENT>(t,agent), nObserved, pObserveIfPresent);
                    assert(observations.back().constraint().isValidSolution(*this));
                }
            }
        }
    }


    double logProb() const {
        StateTrajectory<AGENT> stateTrajectory(*this);
        double logP = 0.0;
        int tEnd = nTimesteps();
        for (int t = 0; t < tEnd; ++t) {
            for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
                double agentOccupation = stateTrajectory[t][agentId];
                if(agentOccupation > 0) {
                    double agentStateLogP = 0.0;
                    std::vector<double> actPMF = AGENT(agentId).timestep(stateTrajectory[t]);
                    double clampedAgentOccupation = 0.0;
                    do {
                        agentStateLogP = 0.0;
                        for (int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                            double actOccupation = (*this)[Event<AGENT>(t, agentId, actId)];
                            if (actOccupation < 0) actOccupation = 0.0;
                            else if (actOccupation > 1.0) actOccupation = 1.0; // Clamp to Fermionic limits
                            if (actOccupation > 0.0) {
                                clampedAgentOccupation += actOccupation;
                                agentStateLogP += actOccupation * log(actPMF[actId]);
                            }
                        }
                        if (agentStateLogP <= -DBL_MAX) actPMF = AGENT(agentId).marginalTimestep();
                    } while(agentStateLogP <= -DBL_MAX);
                    logP += agentStateLogP;
                    if(clampedAgentOccupation > 1.0) logP += lgamma(clampedAgentOccupation + 1.0); // multiply by multinomial
                }
            }
        }
//        assert(fabs(logP - logProbOld()) < 1e-8);
        return logP;
    }


    Trajectory<AGENT> slice(int fromTimestep, int nTimesteps) const {
        Trajectory<AGENT> slice(nTimesteps);
        int beginIndex = Event(fromTimestep,AGENT(0),0);
        int endIndex = Event(fromTimestep+nTimesteps,AGENT(0),0);
        int sliceBeginIndex = Event(0,AGENT(0),0);
        std::copy(begin()+beginIndex, begin()+endIndex, slice.begin()+sliceBeginIndex);
        return slice;
    }


    static double logProb(const std::vector<double> &X) {
        assert(X.size() % (AGENT::domainSize()*AGENT::actDomainSize()) == 0);
        const Trajectory<AGENT> &T = reinterpret_cast<const Trajectory<AGENT> &>(X);
        return T.logProb();
    }

//    friend std::ostream &operator <<(std::ostream &out, const Trajectory<AGENT> &exactEndState) {
//        Trajectory sortedTrajectory(exactEndState);
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

    static std::function<Trajectory<AGENT>()> priorSampler(int nTimesteps, std::function<ModelState<AGENT>()> startStateSampler) {
        return [startStateSampler = std::move(startStateSampler),nTimesteps]() {
            return Trajectory<AGENT>(nTimesteps, startStateSampler);
        };
    }


    static std::vector<double> marginalLogProbsByEvent(int nTimesteps) {
        std::vector<double> marginals(Trajectory<AGENT>::dimension(nTimesteps));
        marginals[0] = 0.0;
        for(int t=0; t<nTimesteps; ++t) {
            for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
                std::vector<double> agentMarginals = AGENT(agentId).marginalTimestep();
                for(int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                    marginals[Event<AGENT>(t,agentId, actId)] = log(agentMarginals[actId]);
                }
            }
        }
        return marginals;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & static_cast<std::vector<double> &>(*this);
    }


};

template<typename AGENT>
typename Trajectory<AGENT>::value_type Trajectory<AGENT>::backwardOccupationNumber(const State<AGENT> &agentState) const {
    assert(agentState.time > 0 && agentState.time <= nTimesteps());
    value_type occupation = 0;
    for(const Event<AGENT> &incomingEvent : State<AGENT>::incomingEventsByState[agentState.agent]) {
        occupation += (*this)[Event(agentState.time-1, incomingEvent.agent(), incomingEvent.act())];
    }
    return occupation;
}

template<typename AGENT>
typename Trajectory<AGENT>::value_type Trajectory<AGENT>::forwardOccupationNumber(const State<AGENT> &agentState) const {
    assert(agentState.time >= 0 && agentState.time < nTimesteps());
    value_type occupation = 0;
    for(int act=0; act<AGENT::actDomainSize(); ++act) {
        occupation += (*this)[Event(agentState.time,agentState.agent,act)];
    }
    return occupation;
}

template<typename AGENT>
typename Trajectory<AGENT>::value_type Trajectory<AGENT>::occupationNumber(const State<AGENT> &agentState) const {
    if(agentState.time == nTimesteps()) return backwardOccupationNumber(agentState);
    return forwardOccupationNumber(agentState);
}


//template<typename AGENT>
//std::function<ModelState<AGENT>()> endStateAdaptor(const std::function<const Trajectory<AGENT> &()> &trajectorySampler) {
//    return [&trajectorySampler]() { return trajectorySampler().endState(); };
//}

#endif //GLPKTEST_TRAJECTORY_H
