//
// Created by daniel on 12/09/22.
//

#ifndef ABMCMC_FIXEDPOPULATIONSTARTSTATE_H
#define ABMCMC_FIXEDPOPULATIONSTARTSTATE_H

#include "ModelState.h"
#include "Trajectory.h"
#include "ConstrainedFactorisedDistribution.h"

template<class AGENT>
class FixedPopulationStartState : public ConstrainedFactorisedDistribution<ModelState<AGENT>> {
public:
    using ConstrainedFactorisedDistribution<ModelState<AGENT>>::constraints;

    explicit FixedPopulationStartState(std::function<int(AGENT)> agentToPartition, std::function<int(int)> partitionToPopulation) {
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            int partitionId = agentToPartition(agentId);
            if(constraints.size() < partitionId+1) constraints.resize(partitionId+1);
            constraints[partitionId].coefficients.insert(agentId,1);
        }
        for(int partitionId=0; partitionId < constraints.size(); ++partitionId) {
            constraints[partitionId].constant = partitionToPopulation(partitionId);
        }
    }


    explicit FixedPopulationStartState(std::initializer_list<int> agentToPartition, std::initializer_list<int> partitionToPopulation):
    FixedPopulationStartState(
            [agentToPartition](AGENT agent) {
                return data(agentToPartition)[agent];
            },
            [partitionToPopulation](int partitionId) {
                return data(partitionToPopulation)[partitionId];
            }) {}


    operator ConstrainedFactorisedDistribution<Trajectory<AGENT>>() {
        return this->toTrajectoryDistribution(0);
    }


    std::function<const ModelState<AGENT> &()> sampler() const {
        return [this, sample = ModelState<AGENT>()]() mutable -> const ModelState<AGENT> & {
            sample = ModelState<AGENT>::zero;
            for(const auto &constraint: constraints) {
                int partitionSize = constraint.coefficients.sparseSize();
                for(int nAgentsPlaced = 0; nAgentsPlaced < constraint.constant; ++nAgentsPlaced) {
                    int chosenPlacement = Random::nextInt(0, partitionSize);
                    int agentId = constraint.coefficients.indices[chosenPlacement];
                    ++sample[agentId];
                }
            }
            return sample;
        };
    }

//    const ModelState<AGENT> &nextSample() const {
//        static thread_local ModelState<AGENT> sample;
//        sample = ModelState<AGENT>::zero;
//        for(const auto &constraint: constraints) {
//            int partitionSize = constraint.coefficients.sparseSize();
//            for(int nAgentsPlaced = 0; nAgentsPlaced < constraint.constant; ++nAgentsPlaced) {
//                int chosenPlacement = Random::nextInt(0, partitionSize);
//                int agentId = constraint.coefficients.indices[chosenPlacement];
//                ++sample[agentId];
//            }
//        }
//        return sample;
//    }


//private:
//    friend class boost::serialization::access;
//
//    template <typename Archive>
//    void serialize(Archive &ar, const unsigned int version) {
//        ar & this->constraints;
//    }

};


#endif //ABMCMC_FIXEDPOPULATIONSTARTSTATE_H
