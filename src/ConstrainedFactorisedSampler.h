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
            ConstrainedFactorisedDistribution<DOMAIN,ELEMENT> &targetDistribution,
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
        std::vector<IndexSet> domainToBasisDependencies;
        std::tie(basisVectors, domainToBasisDependencies) = calculateBasis(targetDistribution.constraints, initialSolution.size());


        // calculate factor by basis-dependency matrix
        std::set<std::pair<int,int>> dependencyMatrix; // entries
        for(const auto &factor: targetDistribution.logFactors) {
            for(int domainIndex: factor.dependencies) {
                for(int basisIndex: domainToBasisDependencies[domainIndex]) {
                    dependencyMatrix.insert(std::pair(factors.size()-1, basisIndex));
                }
            }
        }
        for(auto &dependencyEntry: dependencyMatrix) {
            basisToFactorDependencies[dependencyEntry.second].push_back(dependencyEntry.first);
            factorToBasisDependencies[dependencyEntry.first].push_back(dependencyEntry.second);
        }

        // snap initialSolution to the subspace defined by the basis by recalculating basic vars
        // TODO: Dammit need the TableauNormMinimiser here to retreive F 

        // initialise factors and factor values
        factors.reserve(targetDistribution.logFactors.size());
        currentFactorVal.reserve(targetDistribution.logFactors.size());
        for(const auto &factor: targetDistribution.logFactors) {
            factors.push_back(factor);
            currentFactorVal.push_back(factor(X));
            currentInfeasibility += currentFactorVal.back().second;
        }

        // initialise weights and probabilities
        // TODO: double up the basisDistribution for + and - each basis
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
    static std::pair<std::vector<SparseVec<ELEMENT>>, std::vector<IndexSet>>
    calculateBasis(const EqualityConstraints<ELEMENT> &constraints, int domainDimension) {
        std::pair<std::vector<SparseVec<ELEMENT>>, std::vector<IndexSet>> result;
        std::vector<SparseVec<ELEMENT>> &pmBasisVectors(result.first);
        std::vector<IndexSet> &domainToBasisDependencies(result.second);
        std::vector<SparseVec<ELEMENT>> basisVectors;

        basisVectors = TableauNormMinimiser<ELEMENT>(constraints).getBasisVectors(domainDimension);
        // double-up basis vectors for + and -
        int nBasisVectors = basisVectors.size();
        pmBasisVectors.reserve(nBasisVectors*2);
        for(int i=0; i < basisVectors.size(); ++i) {
            pmBasisVectors.push_back(-basisVectors[i]);
            pmBasisVectors.push_back(std::move(basisVectors[i]));
        }

        // now construct the reverse dependencies by domain index and basis index
        domainToBasisDependencies.resize(domainDimension);
        for(int j=0; j<pmBasisVectors.size(); ++j) {
            for(int row: pmBasisVectors[j].indices) {
                domainToBasisDependencies[row].push_back(j);
            }
        }
        return result;
    }


    const DOMAIN &nextSample() {
        do {
            bool startStateIsFeasible = currentInfeasibility == 0;
            int proposedBasisIndex = basisDistribution(Random::gen);
            double oldSum = basisDistribution.sum();
            performTransition(proposedBasisIndex);
            double acceptance = oldSum / basisDistribution.sum();
            bool proposalIsFeasible = currentInfeasibility == 0;
            bool wasAccepted = true;
            if (Random::nextDouble() > acceptance) {
                wasAccepted = false;
                performTransition(proposedBasisIndex ^ 1); // reject: reverse transition
            }
            stats.addSample(wasAccepted, startStateIsFeasible, proposalIsFeasible);
        } while(currentInfeasibility != 0);
        return X;
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
            currentInfeasibility += currentFactorVal[factorIndex].second - oldLogF.second;
            for(int basisIndex: factorToBasisDependencies[factorIndex]) {
                // update the given factor of the given basis transition
                X += basisVectors[basisIndex];
                double perturbedLogFVal = factors[factorIndex](X).first;
                X -= basisVectors[basisIndex];
                currentLogWeight[basisIndex] += oldLogF.first - 2.0 * currentFactorVal[factorIndex].first + perturbedLogFVal;
                updatedBasisWeights.insert(basisIndex);
            }
        }
        for(int basisIndex: updatedBasisWeights) {
            basisDistribution[basisIndex] = logWeighttoBasisProb(currentLogWeight[basisIndex]);
        }
    }
};


#endif //ABMCMC_CONSTRAINEDFACTORISEDSAMPLER_H
