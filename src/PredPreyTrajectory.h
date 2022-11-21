// An extended trajectory with the addition of variables for the
// count of agents surrounding each gridsquare, suitable for
// predator-prey interaction
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
    PredPreyTrajectory(const SparseVec<int> &sparseInit): ExtendedTrajectory<agent_type,NTIMESTEPS>(dimension, sparseInit) {}

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
//        EqualityConstraints<value_type> constraints = ExtendedTrajectory<agent_type, NTIMESTEPS>::constraints();
        EqualityConstraints<value_type> constraints;
        for (int t = 0; t < NTIMESTEPS; ++t) {
            for (int agentId = 0; agentId < agent_type::domainSize; ++agentId) {
                State<agent_type> state(t, agentId);
//                constraints.push_back(
//                        X(indexOf(State<agent_type>(t, state.agent.upOther()))) +
//                        X(indexOf(State<agent_type>(t, state.agent.downOther()))) +
//                        X(indexOf(State<agent_type>(t, state.agent.leftOther()))) +
//                        X(indexOf(State<agent_type>(t, state.agent.rightOther()))) -
//                        X(surroundingCountIndexOf(state))
//                        == value_type(0)
//                );

                // TODO: need to deal with doubling of birth event when size = 3;
                if(t == 0) {
                    constraints.push_back(
                            X(indexOf(State<agent_type>(t, state.agent.upOther()))) +
                            X(indexOf(State<agent_type>(t, state.agent.downOther()))) +
                            X(indexOf(State<agent_type>(t, state.agent.leftOther()))) +
                            X(indexOf(State<agent_type>(t, state.agent.rightOther()))) -
                            X(surroundingCountIndexOf(state))
                            == value_type(0)
                    );
                } else {
                    SparseVec<value_type> backwardCoeffs;
                    for (const Event<agent_type> &inEdge: State<agent_type>(t,
                                                                            state.agent.upOther()).backwardOccupationDependencies()) {
                        backwardCoeffs.insert(indexOf(inEdge), 1);
                    }
                    for (const Event<agent_type> &inEdge: State<agent_type>(t,
                                                                            state.agent.downOther()).backwardOccupationDependencies()) {
                        backwardCoeffs.insert(indexOf(inEdge), 1);
                    }
                    for (const Event<agent_type> &inEdge: State<agent_type>(t,
                                                                            state.agent.leftOther()).backwardOccupationDependencies()) {
                        backwardCoeffs.insert(indexOf(inEdge), 1);
                    }
                    for (const Event<agent_type> &inEdge: State<agent_type>(t,
                                                                            state.agent.rightOther()).backwardOccupationDependencies()) {
                        backwardCoeffs.insert(indexOf(inEdge), 1);
                    }
                    backwardCoeffs.insert(surroundingCountIndexOf(state), -1);
                    backwardCoeffs.sortAndMerge(); // ensure all entries have unique index
                    constraints.template emplace_back(backwardCoeffs, 0);
                }

            }
        }
        constraints += ExtendedTrajectory<agent_type, NTIMESTEPS>::constraints();
        return constraints;
    }
};
#endif //ABMCMC_PREDPREYTRAJECTORY_H
