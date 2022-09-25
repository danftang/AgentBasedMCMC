//
// Created by daniel on 22/09/22.
//

#ifndef ABMCMC_PREDPREYTRAJECTORY_H
#define ABMCMC_PREDPREYTRAJECTORY_H

#include "ExtendedTrajectory.h"
#include "agents/PredPreyAgent.h"

template<int GRIDSIZE, int NTIMESTEPS>
class PredPreyTrajectory: public ExtendedTrajectory<PredPreyAgent<GRIDSIZE>,NTIMESTEPS> {
public:
    typedef PredPreyAgent<GRIDSIZE> agent_type;
    typedef typename ExtendedTrajectory<agent_type, NTIMESTEPS>::value_type value_type;

    static constexpr size_t dimension =
            ExtendedTrajectory<agent_type, NTIMESTEPS>::dimension +
            NTIMESTEPS * agent_type::domainSize;

    using ExtendedTrajectory<PredPreyAgent<GRIDSIZE>,NTIMESTEPS>::operator [];
    using ExtendedTrajectory<PredPreyAgent<GRIDSIZE>,NTIMESTEPS>::indexOf;


    PredPreyTrajectory() : ExtendedTrajectory<agent_type, NTIMESTEPS>(dimension) {}

    const value_type &surroundingCountOf(const State<agent_type> &state) const {
        return (*this)[surroundingCountIndexOf(state)];
    }

    static int surroundingCountIndexOf(const State<agent_type> &state) {
        return ExtendedTrajectory<agent_type, NTIMESTEPS>::dimension + state.time * agent_type::domainSize +
               state.agent;
    }


    void recalculateDependentVariables(int timestep) {
        ExtendedTrajectory<agent_type, NTIMESTEPS>::recalculateDependentVariables(timestep);
        recalculateSurroundingCounts(timestep);
    }

    void setStartState(const ModelState<PredPreyAgent<GRIDSIZE>> &startState) {
        ExtendedTrajectory<agent_type, NTIMESTEPS>::setStartState(startState);
        recalculateSurroundingCounts(0);
    }

    void recalculateSurroundingCounts(int timestep) {
        assert(timestep >=0 && timestep < NTIMESTEPS);
        for (int agentId = 0; agentId < agent_type::domainSize; ++agentId) {
            State<agent_type> state(timestep, agentId);
            (*this)[surroundingCountIndexOf(state)] =
                    (*this)[State<agent_type>(timestep, state.agent.upOther())] +
                    (*this)[State<agent_type>(timestep, state.agent.downOther())] +
                    (*this)[State<agent_type>(timestep, state.agent.leftOther())] +
                    (*this)[State<agent_type>(timestep, state.agent.rightOther())];
        }
    }


    static EqualityConstraints<value_type> constraints() {
        EqualityConstraints<value_type> constraints = ExtendedTrajectory<agent_type, NTIMESTEPS>::constraints();
        for (int t = 0; t < NTIMESTEPS; ++t) {
            for (int agentId = 0; agentId < agent_type::domainSize; ++agentId) {
                State<agent_type> state(t, agentId);
                constraints.push_back(
                        X(indexOf(State<agent_type>(t, state.agent.upOther()))) +
                        X(indexOf(State<agent_type>(t, state.agent.downOther()))) +
                        X(indexOf(State<agent_type>(t, state.agent.leftOther()))) +
                        X(indexOf(State<agent_type>(t, state.agent.rightOther()))) -
                        X(surroundingCountIndexOf(state))
                        == value_type(0)
                );
            }
        }
        return constraints;
    }
};
#endif //ABMCMC_PREDPREYTRAJECTORY_H
