//
// Created by daniel on 12/09/22.
//

#ifndef ABMCMC_FIXEDPOPULATIONSTARTSTATE_H
#define ABMCMC_FIXEDPOPULATIONSTARTSTATE_H

#include "ModelState.h"
#include "Trajectory.h"
#include "ConstrainedFactorisedDistribution.h"

template<class DOMAIN, class AGENT = typename DOMAIN::agent_type>
class FixedPopulationStartState : public ConstrainedFactorisedDistribution<DOMAIN> {
public:
    using ConstrainedFactorisedDistribution<DOMAIN>::constraints;

    std::function<const ModelState<AGENT> &()> modelStateSampler;

    explicit FixedPopulationStartState(std::function<int(AGENT)> agentToPartition, std::function<int(int)> partitionToPopulation) {
        std::vector<std::vector<int>> partitionToAgents;
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            int partitionId = agentToPartition(agentId);
            if(constraints.size() < partitionId+1) {
                constraints.resize(partitionId+1);
                partitionToAgents.resize(partitionId+1);
            }
            partitionToAgents[partitionId].push_back(agentId);
            State<AGENT> state(0,agentId);
            for(int dependency: DOMAIN::dependencies(state)) {
                constraints[partitionId].coefficients.insert(dependency,1);
            }
        }
        for(int partitionId=0; partitionId < constraints.size(); ++partitionId) {
            constraints[partitionId].constant = partitionToPopulation(partitionId);
        }
        modelStateSampler = sampler(partitionToAgents, partitionToPopulation);
    }


    explicit FixedPopulationStartState(std::initializer_list<int> agentToPartition, std::initializer_list<int> partitionToPopulation):
    FixedPopulationStartState(
            [agentToPartition](AGENT agent) {
                return data(agentToPartition)[agent];
            },
            [partitionToPopulation](int partitionId) {
                return data(partitionToPopulation)[partitionId];
            }) {}


//    operator ConstrainedFactorisedDistribution<Trajectory<AGENT>>() {
//        return this->toTrajectoryDistribution(0);
//    }

    static std::function<const ModelState<AGENT> &()> sampler(std::vector<std::vector<int>> &partitionToAgents, std::function<int(int)> partitionToPopulation) {
        return [partitionToAgents, partitionToPopulation, sample = ModelState<AGENT>()] () mutable -> const ModelState<AGENT> & {
            sample = ModelState<AGENT>::zero;
            for(int partitionId=0; partitionId < partitionToAgents.size(); ++partitionId) {
                int partitionSize = partitionToAgents[partitionId].size();
                int population = partitionToPopulation(partitionId);
                for(int nAgentsPlaced = 0; nAgentsPlaced < population; ++nAgentsPlaced) {
                    int chosenPlacement = Random::nextInt(0, partitionSize);
                    int agentId = partitionToAgents[partitionId][chosenPlacement];
                    ++sample[agentId];
                }
            }
            return sample;
        };
    }

//    std::function<const ModelState<AGENT> &()> sampler() const {
//        return [this, sample = ModelState<AGENT>()]() mutable -> const ModelState<AGENT> & {
//            sample = ModelState<AGENT>::zero;
//            for(const auto &constraint: constraints) {
//                int partitionSize = constraint.coefficients.sparseSize();
//                for(int nAgentsPlaced = 0; nAgentsPlaced < constraint.constant; ++nAgentsPlaced) {
//                    int chosenPlacement = Random::nextInt(0, partitionSize);
//                    int agentId = constraint.coefficients.indices[chosenPlacement];
//                    ++sample[agentId];
//                }
//            }
//            return sample;
//        };
//    }

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
