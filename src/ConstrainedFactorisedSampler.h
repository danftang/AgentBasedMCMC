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
#include "include/debug.h"
#include "include/MutableCategoricalArray.h"
#include "include/SparseVec.h"
#include "include/Random.h"

#include "MCMCStatistics.h"
#include "ConstrainedFactorisedDistribution.h"
#include "TableauNormMinimiser.h"
#include "Basis.h"

template<typename DOMAIN, typename ELEMENT>
class ConstrainedFactorisedSampler {
public:
    typedef std::vector<int> IndexSet;

    struct RowEntry {
        RowEntry(int index): basisIndex(index) { }
        int basisIndex;
        std::pair<double,bool> value;
    };

    struct ColEntry {
        ColEntry(int index, std::pair<double,bool> &val): factorIndex(index), value(val) { }
        int factorIndex;
        std::pair<double,bool> &value;
    };

    std::vector<SparseVec<ELEMENT>> basisVectors;       // the basis vectors by basis index and domain index
    MutableCategoricalArray         basisDistribution;  // probability of proposing to update the j'th basis
    std::vector<double>             currentWeight;   // current log weight of the j'th transition

    std::vector<std::function<std::pair<double,bool>(const DOMAIN &)>>  factors;
    std::vector<std::pair<double,bool>>                         currentFactorVal;       // current log value of the i'th factor

    // dependency matrix, each col is a basis, each row is a factor
    std::vector<std::vector<RowEntry>>       factorToBasisPerturbedValue; // set of basis vectors that a given factor is dependent on
    std::vector<std::vector<ColEntry>>     basisToFactorPerturbedValue; // set of factors that are dependent on a given basis

    DOMAIN                      X;      // the current point on the lattice
    int                         currentInfeasibility; // current number of factors that are in widened support
    MCMCStatistics              stats;

public:
    ConstrainedFactorisedSampler(const ConstrainedFactorisedDistribution<DOMAIN, ELEMENT> &targetDistribution):
            ConstrainedFactorisedSampler(targetDistribution, Basis(targetDistribution)) { }

    ConstrainedFactorisedSampler(const ConstrainedFactorisedDistribution<DOMAIN, ELEMENT> &targetDistribution, const Basis<DOMAIN> &basis):
            ConstrainedFactorisedSampler(targetDistribution, basis.basisVectors, basis.origin) { }

    ConstrainedFactorisedSampler(
            const ConstrainedFactorisedDistribution<DOMAIN, ELEMENT> &targetDistribution,
            const std::vector<SparseVec<ELEMENT>> &basisVecs,
            const DOMAIN &initialValidSample): X(initialValidSample) {

        assert(targetDistribution.constraints.isValidSolution(X));

        // double-up unique basis vectors for + and -
        basisVectors.reserve(basisVecs.size()*2);
        for(int i=0; i < basisVecs.size(); ++i) {
            basisVectors.push_back(-(basisVecs[i]));
            basisVectors.push_back( basisVecs[i]);
        }

        initDependencies(targetDistribution, basisVectors);

        // initialise factors and factor values
        factors.reserve(targetDistribution.factors.size());
        currentFactorVal.reserve(targetDistribution.factors.size());
        currentInfeasibility = 0;
        for(const auto &factor: targetDistribution.factors) {
            factors.push_back(factor);
            currentFactorVal.push_back(factor(X));
            currentInfeasibility += !currentFactorVal.back().second;
        }

        // initialise weights, ratios and probabilities
        currentWeight.resize(basisVectors.size(), 0.0);
        for(int basisIndex=0; basisIndex < basisVectors.size(); ++basisIndex) {
            X += basisVectors[basisIndex];
            for(ColEntry &colEntry : basisToFactorPerturbedValue[basisIndex]) {
                colEntry.value = factors[colEntry.factorIndex](X);
                currentWeight[basisIndex] += colEntry.value.first - currentFactorVal[colEntry.factorIndex].first;
            }
            X -= basisVectors[basisIndex];
            basisDistribution.push_back(logWeighttoBasisProb(currentWeight[basisIndex]));
        }
        debug(sanityCheck());
    }

    const DOMAIN &operator()() { return nextSample(); }

    const DOMAIN &nextSample() {
        debug(int attempts = 0);
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
//        double PX0 = currentProb();
        for(int basisId=0; basisId < basisVectors.size(); ++basisId) {
            X += basisVectors[basisId];
            double logRatioOfProbs = currentLogProb() - logProbX0;
            assert(fabs(currentWeight[basisId] - logRatioOfProbs) < 1e-8);
            assert(fabs(logWeighttoBasisProb(logRatioOfProbs) - basisDistribution[basisId]) < 1e-8);
//            double ratioOfProbs = currentProb()/PX0;
//            assert(fabs(currentWeight[basisId] - ratioOfProbs) < 1e-8);
//            assert(fabs(weighttoBasisProb(ratioOfProbs) - basisDistribution[basisId]) < 1e-8);
            X -= basisVectors[basisId];
        }
    }


protected:

    // Turns the given constraints into a minimal basis that spans the
    // subspace described by the constraints.
    // Returns a pair whose first element is the basis vectors
    // and whose second element is the dependency mapping from domain element to
    // the set of basis vectors that depend on that element, where the
    // basis vectors are given an index by the ordering of the first element.
    void initDependencies(const ConstrainedFactorisedDistribution<DOMAIN,ELEMENT> &targetDistribution,
                          const std::vector<SparseVec<ELEMENT>> &basisVecs)
    {

        // construct the dependencies from domain index to basis index
        std::vector<IndexSet> domainToBasisDependencies(DOMAIN::dimension);
        int nBasisEntries = 0;
        for(int j=0; j<basisVecs.size(); ++j) {
            for(int row: basisVecs[j].indices) {
                domainToBasisDependencies[row].push_back(j);
                ++nBasisEntries;
            }
        }

        // calculate factor by basis-dependency matrix
        std::set<std::pair<int,int>> dependencyMatrix; // entries
        for(int factorId = 0; factorId<targetDistribution.factors.size(); ++factorId) {
            for(int domainIndex: targetDistribution.factors[factorId].dependencies) {
                for(int basisIndex: domainToBasisDependencies[domainIndex]) {
                    dependencyMatrix.insert(std::pair(factorId, basisIndex));
                }
            }
        }
        std::cout << "mean basis size = " << nBasisEntries * 1.0 / basisVecs.size() << std::endl;
        std::cout << "mean dependency col size = " << dependencyMatrix.size() * 1.0 / basisVecs.size() << std::endl;
        std::cout << "mean dependency row size = " << dependencyMatrix.size() * 1.0/ targetDistribution.factors.size() << std::endl;
        std::cout << "mean updates per transition = " << pow(dependencyMatrix.size(),2)*1.0/(basisVecs.size()*targetDistribution.factors.size())  << std::endl;

        // now transfer to forward and backward dependencies
        // set up structure, values will be added later
        factorToBasisPerturbedValue.resize(targetDistribution.factors.size());
        for(auto &dependencyEntry: dependencyMatrix) {
            factorToBasisPerturbedValue[dependencyEntry.first].emplace_back(dependencyEntry.second);
        }
        basisToFactorPerturbedValue.resize(basisVecs.size());
        for(int factorId = 0; factorId < factorToBasisPerturbedValue.size(); ++factorId) {
            for(RowEntry &entry : factorToBasisPerturbedValue[factorId]) {
                basisToFactorPerturbedValue[entry.basisIndex].emplace_back(factorId, entry.value);
            }
        }
    }

    static double logWeighttoBasisProb(double w) { return w<0.0?exp(w):1.0; } // convert change in energy of a column to probability of proposal


    // Given a transition along basis j,
    // Update each basis, k, for which there exists a factor that depends on both j and k
    void performTransition(int transitionBasis) {
        std::set<int> updatedBasisWeights;
        X += basisVectors[transitionBasis];
        for(ColEntry &colEntry: basisToFactorPerturbedValue[transitionBasis]) {
            int factorIndex = colEntry.factorIndex;
            double newFVal = colEntry.value.first;
            double oldFVal = currentFactorVal[factorIndex].first;
            for(RowEntry &rowEntry: factorToBasisPerturbedValue[factorIndex]) {
                int basisIndex = rowEntry.basisIndex;
                std::pair<double,bool> &factorVal = rowEntry.value;
                if(basisIndex == transitionBasis) {
                    currentInfeasibility += currentFactorVal[factorIndex].second - factorVal.second;
                    currentFactorVal[factorIndex] = factorVal;
                }
                currentWeight[basisIndex] += oldFVal - factorVal.first; // remove old ratio
                // update the given factor of the given basis transition
                X += basisVectors[basisIndex];
                factorVal = factors[factorIndex](X);
                X -= basisVectors[basisIndex];
                currentWeight[basisIndex] += factorVal.first - newFVal; // insert new ratio
                updatedBasisWeights.insert(basisIndex);
            }
        }
        for(int basisIndex: updatedBasisWeights) {
            basisDistribution[basisIndex] = logWeighttoBasisProb(currentWeight[basisIndex]);
        }
    }
};


#endif //ABMCMC_CONSTRAINEDFACTORISEDSAMPLER_H
