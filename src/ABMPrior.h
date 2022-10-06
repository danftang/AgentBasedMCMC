// Represents the prior distribution over an Agent based model
//
// Created by daniel on 01/04/2022.
//

#ifndef ABMCMC_ABMPRIOR_H
#define ABMCMC_ABMPRIOR_H

#include "ABM.h"
#include "ConstrainedFactorisedDistribution.h"

// DOMAIN should be a domain over ABM trajectories that inplements
// index operators over Events and States
template<class TRAJECTORY, class STARTSTATE, class AGENT = typename TRAJECTORY::agent_type>
class ABMPrior: public ConstrainedFactorisedDistribution<TRAJECTORY> {
public:
//    std::function<const ModelState<AGENT> &()> startStateSampler;
    STARTSTATE startState;

    ABMPrior()=default;

    ABMPrior(STARTSTATE startState): startState(std::move(startState)) {
        init();
    }

protected:
    void init() {
        this->addConstraints(TRAJECTORY::constraints());
        addMultinomials();
        (*this) *= startState.template toDomain<TRAJECTORY>();
    }

public:

    void addMultinomials() {
        for(int time = 0; time < TRAJECTORY::nTimesteps; ++time) {
            for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
                State<AGENT> state(time, AGENT(agentId));

                for (int actId = 0; actId < AGENT::actDomainSize; ++actId) {
                    Event<typename TRAJECTORY::agent_type> event(time, agentId, actId);
                    std::vector<int> dependencies = TRAJECTORY::agent_type::template eventProbDependencies<TRAJECTORY>(event);
                    dependencies.push_back(TRAJECTORY::indexOf(event));
                    this->addFactor(
                            [event, kappa = startState.kappa](const TRAJECTORY &trajectory) {
                                return widenedEventFactor(event, trajectory, kappa);
                            },
                            dependencies
                    );
                }

                this->addFactor( // factorial of state occupation
                        [state](const TRAJECTORY &trajectory) {
                            auto occupation = trajectory[state];
//                            return (occupation < 0) ? std::pair(ABM::kappa*occupation,false) : std::pair(lgamma(occupation + 1), true); // TODO: Test!!
//                            return std::pair(lgamma(std::abs(trajectory[state]) + 1), true); // TODO: test
                            return std::pair(lgamma(std::max(trajectory[state], 0) + 1), true);
                        },
                        {TRAJECTORY::indexOf(state)}
                );
            }
        }
    }



    static std::pair<double, bool> widenedEventFactor(const Event<AGENT> &event, const TRAJECTORY &trajectory, double kappa) {
        int actOccupation = trajectory[event];
        if(actOccupation == 0) return std::pair(0.0, true);
        if (actOccupation < 0) { // negative occupation widening
//            return std::pair(-actOccupation * log(1.0/AGENT::actDomainSize) - lgamma(1-actOccupation) + ABM::kappa * actOccupation, false); // pair production and decay
//            return std::pair(ABM::kappa * actOccupation - lgamma(1-actOccupation), false);
            return std::pair(kappa * actOccupation, false);
        }

//        auto [logpAct, isFeasible] = AGENT::widenedLogEventProb(event, trajectory); // TODO: test!!!
//        if(actOccupation == 1) return std::pair(logpAct, isFeasible); // just to save doing the logGamma
//        if(actOccupation > 1) { // Fermionic occupation widening
//            return std::pair(logpAct + ABM::kappa * (1-actOccupation), false); // TODO: test!!!
//        }
//        return std::pair(actOccupation * logpAct - lgamma(actOccupation + 1), isFeasible);

        double logpAct = AGENT::logEventProb(event, trajectory);
        if (logpAct == -INFINITY) { // impossible act widening
//            return std::pair(event.agent().logMarginalTimestep(event.act())* actOccupation - lgamma(actOccupation + 1) - kappa, false);
            return std::pair(-kappa*actOccupation - lgamma(actOccupation + 1), false); // TODO: test!!!
        }
        if(actOccupation == 1) return std::pair(logpAct, true); // just to save doing the logGamma
        return std::pair(actOccupation * logpAct - lgamma(actOccupation + 1), true);
    }



    // fills only the neighbours of an agent-state, leaving the other elements undefined
//    static const ModelState<AGENT> &neighbourModelState(const TRAJECTORY &trajectory, const State<AGENT> &agentState) {
//        static thread_local ModelState<AGENT> temporaryModelState;
//        for(AGENT neighbour: agentState.agent.eventProbDependencies())
//            temporaryModelState[neighbour] = trajectory[State<AGENT>(agentState.time, neighbour)];
//        return temporaryModelState;
//    }

    const TRAJECTORY &operator ()() const { return nextSample(); }

    const TRAJECTORY &nextSample() const {
        static thread_local TRAJECTORY sample;
        sample.setStartState(startState.nextSample());
        for (int t = 0; t < TRAJECTORY::nTimesteps; ++t) {
            if (t > 0) sample.recalculateDependentVariables(t);
            for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
                AGENT agent(agentId);
                State<AGENT> state(t, agent);
                int nAgents = sample[state];
                for (int actId = 0; actId < AGENT::actDomainSize; ++actId) {
                    sample[TRAJECTORY::indexOf(Event<AGENT>(t, agent, actId))] = 0; // modification must be via int index
                }
                std::vector<double> actPMF = agentActPMF(state, sample);
                for (int a = 0; a < nAgents; ++a) {
                    int nextAct = Random::nextIntFromDiscrete(actPMF);
                    sample[TRAJECTORY::indexOf(Event<AGENT>(t, agent, nextAct))] += 1;
                }
            }
        }
        return sample;
    }


    static std::vector<double> agentActPMF(const State<AGENT> &state, const TRAJECTORY &trajectory) {
        std::vector<double> pmf;
        pmf.reserve(AGENT::actDomainSize);
        for(int actId=0; actId < AGENT::actDomainSize; ++actId) {
            Event<AGENT> event(state.time, state.agent, actId);
            pmf.push_back(exp(AGENT::logEventProb(event, trajectory)));
        }
        return pmf;
    }

    friend std::ostream &operator <<(std::ostream &out, const ABMPrior<TRAJECTORY,STARTSTATE> &prior) {
        out << "Start state: " << prior.startState << std::endl;
        out << "Prior: " << static_cast<ConstrainedFactorisedDistribution<TRAJECTORY>>(prior);
        return out;
    }

private:

    friend class boost::serialization::access;

    template <typename Archive>
    void load(Archive &ar, const unsigned int version) {
        ar >> startState;
        init();
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        ar << startState;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();


    // Turns an agent timestep, pi(a, \Psi), unction into a function
    // \prod_a \pi(a, \Psi, \psi)^T^t_{\psi a}/T^t_{\psi a}!
    // widened to decay (roughly) exponentially over negative occupations and zero probability actions
//    static std::pair<double, bool> widenedAgentMultinomial(const State<AGENT> &state, const DOMAIN &trajectory, const bool includePhiFactorial) {
//        ABM::occupation_type stateOccupation(0);
//        double logP(0.0);
//        bool exactValue = true;
//        // std::vector<double> actPMF = state.agent.timestep(
////        const ModelState<AGENT> &neighbours = neighbourModelState(trajectory, state);
//        for (int act = 0; act < AGENT::actDomainSize; ++act) {
//            Event<AGENT> event(state.time, state.agent, act);
//            int actOccupation = trajectory[event];
//            if (actOccupation != 0) {
//                double logpAct;
//                if (actOccupation < 0) { // negative occupation widening
//                    logpAct = -ABM::kappa;
//                    actOccupation = -actOccupation;
//                    exactValue = false;
//                } else {
//                    logpAct = AGENT::logEventProb(event, trajectory);
//                    if (logpAct == -INFINITY) {
//                        logpAct = -ABM::kappa; // impossible act widening
//                        exactValue = false;
//                    }
//                }
//                logP += actOccupation * logpAct - lgamma(actOccupation + 1);
//                stateOccupation += actOccupation;
//            }
//        }
//        if(includePhiFactorial) logP += lgamma(std::max(stateOccupation, 0) + 1); // Phi factorial
//        return std::pair(logP, exactValue);
//    }

    //    std::function<const DOMAIN &()> sampler() const {
//        return [sample = DOMAIN(), startStateSampler = startStateSampler]() mutable -> const DOMAIN & {
//            sample.setStartState(startStateSampler());
//            for (int t = 0; t < DOMAIN::nTimesteps; ++t) {
//                if(t>0) sample.recalculateDependentVariables(t);
//                for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
//                    AGENT agent(agentId);
//                    State<AGENT> state(t, agent);
//                    int nAgents = sample[state];
//                    for (int actId = 0; actId < AGENT::actDomainSize; ++actId) {
//                        sample[DOMAIN::indexOf(Event<AGENT>(t, agent, actId))] = 0; // modification must be via int index
//                    }
//                    std::vector<double> actPMF = agentActPMF(state, sample);
//                    for (int a = 0; a < nAgents; ++a) {
//                        int nextAct = Random::nextIntFromDiscrete(actPMF);
//                        sample[DOMAIN::indexOf(Event<AGENT>(t, agent, nextAct))] += 1;
//                    }
//                }
//            }
//            return sample;
//        };
//    }


//    friend std::ostream &operator <<(std::ostream &out, const Prior &prior) {
//        out << "Prior timesteps = " << prior.nTimesteps << "\nPrior ModelState =\n" << prior.startState << std::endl;
//        out << "Prior Trajectory constraints =\n" << prior.constraints;
//        return out;
//    }

};

template<class TRAJECTORY, class STARTSTATE>
ABMPrior<TRAJECTORY,STARTSTATE> makeABMPrior(STARTSTATE startstate) { return ABMPrior<TRAJECTORY,STARTSTATE>(startstate); }

#endif //ABMCMC_ABMPRIOR_H
