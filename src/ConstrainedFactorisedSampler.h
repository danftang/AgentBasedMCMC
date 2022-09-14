// Sampler for ConstrainedFactoriseedDistribution
// Represents a Markov process whose current state is X, and whose
// transitions are a set of basis vectors which are added to X
//
// Each state has an un-normalised probability and each
// transition has a probability
//
// IDEA: Suppose current state is X and we're considering
// a transition to X'. If P(X') > P(X) then transition
// with probability 1/N where N is the number of
// basis vectors. If P(X') < P(X) then transition
// with probability P(X')/(NP(X)) the resultant
// sum of transition probabilities are guaranteed to
// be no greater than 1, so fill any remainder with
// the loop transition. We then only need to know the
// sum of transitions of the current state and we don't
// need to deal with rejection.
//
// Created by daniel on 07/09/22.
//

#ifndef ABMCMC_CONSTRAINEDFACTORISEDSAMPLER_H
#define ABMCMC_CONSTRAINEDFACTORISEDSAMPLER_H

#include <vector>
#include <map>
#include <functional>
#include "MutableCategoricalArray.h"
#include "SparseVec.h"
#include "MCMCStatistics.h"
#include "Random.h"
#include "ConstrainedFactorisedDistribution.h"
#include "TableauNormMinimiser.h"

template<typename DOMAIN, typename ELEMENT = typename subscript_operator_traits<DOMAIN>::base_type>
class ConstrainedFactorisedSampler {
public:
    typedef std::vector<int> IndexSet;

    std::vector<SparseVec<ELEMENT>> basisVectors;       // the basis vectors by basis index and domain index
    MutableCategoricalArray         basisDistribution;  // probability of proposing to update the j'th basis
    std::vector<double>             currentLogWeight;   // current log weight of the j'th transition

    std::vector<std::function<std::pair<double,bool>(const DOMAIN &)>>  factors;
    std::vector<std::pair<double,bool>>                         currentFactorVal;       // current log value of the i'th factor

    // dependency matrix, each col is a basis, each row is a factor
    std::vector<IndexSet>       factorToBasisDependencies; // set of basis vectors that a given factor is dependent on
    std::vector<IndexSet>       basisToFactorDependencies; // set of factors that are dependent on a given basis

    DOMAIN                      X;      // the current point on the lattice by row.
    int                         currentInfeasibility; // current number of factors that are in widened support
    MCMCStatistics              stats;

    ConstrainedFactorisedSampler<DOMAIN,ELEMENT>(
            const ConstrainedFactorisedDistribution<DOMAIN,ELEMENT> &targetDistribution,
            DOMAIN initialSolution
            ):
            basisDistribution(2*(initialSolution.size() - targetDistribution.constraints.size())),
            currentLogWeight(basisDistribution.size()),
            factorToBasisDependencies(targetDistribution.logFactors.size()),
            basisToFactorDependencies(basisDistribution.size()),
            X(initialSolution),
            currentInfeasibility(0)
    {

        // calculate basis
        TableauNormMinimiser<ELEMENT> constraintTableau(targetDistribution.constraints);
        std::vector<SparseVec<ELEMENT>> uniqueBasisVectors = constraintTableau.getBasisVectors(initialSolution.size());
        // double-up unique basis vectors for + and -
        basisVectors.reserve(uniqueBasisVectors.size()*2);
        for(int i=0; i < uniqueBasisVectors.size(); ++i) {
            basisVectors.push_back(-uniqueBasisVectors[i]);
            basisVectors.push_back(std::move(uniqueBasisVectors[i]));
        }
        constraintTableau.snapToSubspace(X);
        assert(targetDistribution.constraints.isValidSolution(X));
        initDependencies(targetDistribution);

        // Test the basis
//        for(int i=0; i<basisVectors.size(); ++i) {
//            X += basisVectors[i];
//            std::cout << "Testing basis vector " << basisVectors[i] << " X = " << X << std::endl;
//            assert(targetDistribution.constraints.isValidSolution(X));
//            X -= basisVectors[i];
//        }

        // initialise factors and factor values
        factors.reserve(targetDistribution.logFactors.size());
        currentFactorVal.reserve(targetDistribution.logFactors.size());
        for(const auto &factor: targetDistribution.logFactors) {
            factors.push_back(factor);
            currentFactorVal.push_back(factor(X));
            currentInfeasibility += !currentFactorVal.back().second;
        }

        // initialise weights and probabilities
        for(int j=0; j<basisVectors.size(); ++j) {
            currentLogWeight[j] = 0.0;
            X += basisVectors[j];
            for(int i: basisToFactorDependencies[j]) {
                currentLogWeight[j] += factors[i](X).first - currentFactorVal[i].first;
            }
            X -= basisVectors[j];
            basisDistribution[j] = logWeighttoBasisProb(currentLogWeight[j]);
        }
    }

    // Turns the given constraints into a minimal basis that spans the
    // subspace described by the constraints.
    // Returns a pair whose first element is the basis vectors
    // and whose second element is the dependency mapping from domain element to
    // the set of basis vectors that depend on that element, where the
    // basis vectors are given an index by the ordering of the first element.
    void initDependencies(const ConstrainedFactorisedDistribution<DOMAIN,ELEMENT> &targetDistribution) {
        std::vector<IndexSet> domainToBasisDependencies(X.size());

        // construct the dependencies from domain index to basis index
        for(int j=0; j<basisVectors.size(); ++j) {
            for(int row: basisVectors[j].indices) {
                domainToBasisDependencies[row].push_back(j);
            }
        }

        // calculate factor by basis-dependency matrix
        std::set<std::pair<int,int>> dependencyMatrix; // entries
        for(int factorId = 0; factorId<targetDistribution.logFactors.size(); ++factorId) {
            for(int domainIndex: targetDistribution.logFactors[factorId].dependencies) {
                for(int basisIndex: domainToBasisDependencies[domainIndex]) {
                    dependencyMatrix.insert(std::pair(factorId, basisIndex));
                }
            }
        }

        // now transfer to forward and backward maps
        for(auto &dependencyEntry: dependencyMatrix) {
            basisToFactorDependencies[dependencyEntry.second].push_back(dependencyEntry.first);
            factorToBasisDependencies[dependencyEntry.first].push_back(dependencyEntry.second);
        }
    }


    const DOMAIN &nextSample() {
        int attempts = 0;
        do {
//            std::cout << "Infeasibility = " << currentInfeasibility << std::endl;
            bool startStateIsFeasible = currentInfeasibility == 0;
            int proposedBasisIndex = basisDistribution(Random::gen);
            double oldSum = basisDistribution.sum();
            performTransition(proposedBasisIndex);
            double acceptance = oldSum / basisDistribution.sum();
            bool proposalIsFeasible = currentInfeasibility == 0;
            bool wasAccepted = true;
            if (Random::nextDouble() > acceptance) {
//                std::cout << "Rejecting" << std::endl;
                wasAccepted = false;
                performTransition(proposedBasisIndex ^ 1); // reject: reverse transition
            } else {
//                std::cout << "Accepting" << std::endl;
            }
            stats.addSample(wasAccepted, startStateIsFeasible, proposalIsFeasible);
            debug(if(++attempts%10000 == 0) std::cout << "Stuck with infeasibility = " << currentInfeasibility << std::endl;);
            debug(sanityCheck());
        } while(currentInfeasibility != 0);
        return X;
    }

    double currentLogProb() {
        double logProb = 0.0;
        for(const auto &factor: factors) logProb += factor(X).first;
        return logProb;
    }


    void sanityCheck() {
        // check weights
        double logProbX0 = currentLogProb();
        for(int basisId=0; basisId < basisVectors.size(); ++basisId) {
            X += basisVectors[basisId];
            double logRatioOfProbs = currentLogProb() - logProbX0;
            assert(fabs(currentLogWeight[basisId] - logRatioOfProbs) < 1e-8);
            assert(fabs(logWeighttoBasisProb(logRatioOfProbs) - basisDistribution[basisId]) < 1e-8);
            X -= basisVectors[basisId];
        }
    }

protected:
    static double logWeighttoBasisProb(double w) { return w<0.0?exp(w):1.0; } // convert change in energy of a column to probability of proposal


    // Given a transition along basis j,
    // Update each basis, k, for which there exists a factor that depends on both j and k
    void performTransition(int perturbedBasis) {
        std::set<int> updatedBasisWeights;
        X += basisVectors[perturbedBasis];
        for(int factorIndex: basisToFactorDependencies[perturbedBasis]) {
            std::pair<double,bool> oldLogF = currentFactorVal[factorIndex];
            currentFactorVal[factorIndex] = factors[factorIndex](X);
            currentInfeasibility += oldLogF.second - currentFactorVal[factorIndex].second;
            for(int basisIndex: factorToBasisDependencies[factorIndex]) {
                // update the given factor of the given basis transition
                X += basisVectors[basisIndex];   // TODO: is it worth storing these?... or nicer to have add/remove(basis, factor)?
                X -= basisVectors[perturbedBasis];
                double oldPerturbedLogFVal = factors[factorIndex](X).first;
                X += basisVectors[perturbedBasis];
                double perturbedLogFVal = factors[factorIndex](X).first;
                X -= basisVectors[basisIndex];
                currentLogWeight[basisIndex] += oldLogF.first - oldPerturbedLogFVal + perturbedLogFVal - currentFactorVal[factorIndex].first;
                updatedBasisWeights.insert(basisIndex);
            }
        }
        for(int basisIndex: updatedBasisWeights) {
            basisDistribution[basisIndex] = logWeighttoBasisProb(currentLogWeight[basisIndex]);
        }
    }
};


#endif //ABMCMC_CONSTRAINEDFACTORISEDSAMPLER_H
