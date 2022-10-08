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

#ifndef ABMCMC_FACTORISEDDISTRIBUTIONSAMPLER_H
#define ABMCMC_FACTORISEDDISTRIBUTIONSAMPLER_H

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
#include "FactorisedDistributionOnBasis.h"
#include "ABMPosterior.h"



template<typename DOMAIN, typename ELEMENT = typename DOMAIN::value_type>
class FactorisedDistributionSampler {
public:

    // USER is a function object that takes a producer of samples and returns a result
//    template<class USER, class RESULT = typename SamplingProblem<ConstrainedFactorisedDistribution<DOMAIN,ELEMENT>,USER>::result_type>
//    static RESULT solve(const SamplingProblem<ConstrainedFactorisedDistribution<DOMAIN,ELEMENT>,USER> &problem) {
//        FactorisedDistributionSampler sampler(problem.distribution);
//        return problem.user(sampler);
//    }



public:
    typedef std::vector<int> IndexSet;

    struct RowEntry {
        RowEntry(int index): basisIndex(index) { }
        int basisIndex;
        std::pair<double,bool> value; // value of this row's factor after unit perturbation of this column's basis vector
    };

    struct ColEntry {
        ColEntry(int index, std::pair<double,bool> &val): factorIndex(index), value(val) { }
        int factorIndex;
        std::pair<double,bool> &value; // reference to value stored in corresponding RowEntry
    };

    std::vector<SparseVec<ELEMENT>> basisVectors;       // the basis vectors by basis index and domain index (doubled for +-)
    std::vector<std::function<std::pair<double,bool>(const DOMAIN &)>>  factors;
    // dependency matrix, each col is a basis, each row is a factor
    std::vector<std::vector<RowEntry>>       factorToBasisPerturbedValue; // set of basis vectors that a given factor is dependent on
    std::vector<std::vector<ColEntry>>     basisToFactorPerturbedValue; // set of factors that are dependent on a given basis

    MutableCategoricalArray         basisDistribution;  // probability of proposing to update the j'th basis
    std::vector<double>             currentWeight;   // current log weight of the j'th transition
    std::vector<std::pair<double,bool>>  currentFactorVal;  // current log value of the i'th factor

    DOMAIN                      X;      // the current Markov state
    int                         currentInfeasibility; // current number of factors that are in widened support
    MCMCStatistics              stats;  // statistics on all samples since construction

public:
    FactorisedDistributionSampler(FactorisedDistributionSampler &&other)=delete;
    FactorisedDistributionSampler(const FactorisedDistributionSampler &other)=delete;
//    basisVectors(other.basisVectors),
//    factors(other.factors),
//    X(other.X) {
//        init(factors, basisVectors, X);
//    }

    FactorisedDistributionSampler(const ConstrainedFactorisedDistribution<DOMAIN, ELEMENT> &targetDistribution):
            FactorisedDistributionSampler(Basis(targetDistribution), targetDistribution.factors) { }
//
////    FactorisedDistributionSampler(const ConstrainedFactorisedDistribution<DOMAIN, ELEMENT> &targetDistribution, const Basis<DOMAIN> &basis):
////            FactorisedDistributionSampler(targetDistribution, basis.basisVectors, basis.origin) { }
//
//    FactorisedDistributionSampler(const FactorisedDistributionOnBasis<DOMAIN, ELEMENT> &targetDistribution):
//            FactorisedDistributionSampler(targetDistribution, targetDistribution.basis.origin) { }
//
//    FactorisedDistributionSampler(const ConstrainedFactorisedDistribution<DOMAIN, ELEMENT> &targetDistribution, const DOMAIN &initialValidSample):
//            FactorisedDistributionSampler(FactorisedDistributionOnBasis(targetDistribution), initialValidSample) { }

    template<class STARTSTATE>
    FactorisedDistributionSampler(const ABMPosterior<DOMAIN,STARTSTATE> &targetDistribution):
            FactorisedDistributionSampler(targetDistribution.basis, targetDistribution.factors) { }


    FactorisedDistributionSampler(
            const Basis<DOMAIN> &basis,
            const std::vector<SparseFunction<std::pair<double,bool>,const DOMAIN &>> &distFactors): X(basis.origin) {

        Random::gen.seed(Random::nextRandomSeed());
        for(int i=0; i<basis.basisVectors.size(); ++i) {
            if(Random::nextBool(0.005)) X += basis.basisVectors[i];
        }

//        assert(targetDistribution.constraints.isValidSolution(X));

        // double-up unique basis vectors for + and -
        basisVectors.reserve(basis.basisVectors.size()*2);
        for(int i=0; i < basis.basisVectors.size(); ++i) {
            basisVectors.push_back(-(basis.basisVectors[i]));
            basisVectors.push_back( basis.basisVectors[i]);
        }

        factors.reserve(distFactors.size());
        for(const auto &factor: distFactors) {
            factors.push_back(factor);
        }

        init(distFactors, basisVectors, X);

//        // initialise factors and factor values
//        factors.reserve(distFactors.size());
//        currentFactorVal.reserve(distFactors.size());
//        currentInfeasibility = 0;
//        for(const auto &factor: distFactors) {
//            factors.push_back(factor);
//            currentFactorVal.push_back(factor(X)); // Assumes all factors with single basis dependency are feasible
//            assert((factorToBasisPerturbedValue[factors.size()-1].size() != 2) || (currentFactorVal.back().second == true));
//            currentInfeasibility += !currentFactorVal.back().second;
//        }
//
//        // initialise weights, ratios and probabilities
//        currentWeight.resize(basisVectors.size(), 0.0);
//        for(int basisIndex=0; basisIndex < basisVectors.size(); ++basisIndex) {
//            X += basisVectors[basisIndex];
//            for(ColEntry &colEntry : basisToFactorPerturbedValue[basisIndex]) {
//                colEntry.value = factors[colEntry.factorIndex](X);
////                std::cout << factorToBasisPerturbedValue[colEntry.factorIndex].size() << std::endl;
//                if(factorToBasisPerturbedValue[colEntry.factorIndex].size() == 2 && colEntry.value.second == false) { // if this factor depends only on this basis, stay within feasible region
//                    colEntry.value.first -= 24.0; // TODO: test!!!!
//                }
//                currentWeight[basisIndex] += colEntry.value.first - currentFactorVal[colEntry.factorIndex].first;
//            }
//            X -= basisVectors[basisIndex];
//            basisDistribution.push_back(logWeighttoBasisProb(currentWeight[basisIndex]));
//        }
//        debug(sanityCheck());
        findInitialFeasibleSolution();
    }

    // X should start valid, but not necessarily feasible
    void findInitialFeasibleSolution() {
        int nTransitions = 0;
        while(currentInfeasibility != 0) {
//            std::cout << "current infeasibility " << currentInfeasibility << std::endl;
            performTransition(basisDistribution(Random::gen));
            ++nTransitions;
        }
        std::cout << "Found initial feasible solution in " << nTransitions << " transitions" << std::endl;
    }

    const DOMAIN &operator()() { return nextSample(); }

    const DOMAIN &nextSample() {
        debug(int attempts = 0);
        do {
//            std::cout << "Infeasibility = " << currentInfeasibility << std::endl;
            bool startStateIsFeasible = (currentInfeasibility == 0);
            int proposedBasisIndex = basisDistribution(Random::gen);
            double oldSum = basisDistribution.sum();
            performTransition(proposedBasisIndex);
            bool proposalIsFeasible = (currentInfeasibility == 0);
            double acceptance = oldSum / basisDistribution.sum();
//            double acceptance = (oldSum*(proposalIsFeasible?1.0:0.9)) / (basisDistribution.sum()*(startStateIsFeasible?1.0:0.9)); // TODO: Test!!!
            bool wasAccepted = true;
            if (Random::nextDouble() > acceptance) {
//                std::cout << "Rejecting" << std::endl;
                wasAccepted = false;
                performTransition(proposedBasisIndex ^ 1); // reject: reverse transition
            } else {
//                std::cout << "Accepting. Infeasibility = " << currentInfeasibility << std::endl;
            }
            stats.addSample(wasAccepted, startStateIsFeasible, proposalIsFeasible);
            debug(if(++attempts%10000 == 0) std::cout << attempts << ": stuck with infeasibility = " << currentInfeasibility << std::endl;);
            debug(sanityCheck());
        } while(currentInfeasibility != 0);
        return X;
    }

    double currentLogProb() {
        double logProb = 0.0;
        for(int i=0; i<factors.size(); ++i) {
            auto factorVal = factors[i](X);
            logProb += factorVal.first;
            if(factorToBasisPerturbedValue[i].size() == 2 && factorVal.second == false) {
                logProb -= 24.0; // TODO: Test!!!!
            }
        }
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
        int infeasibility = 0;
        for(int j=0; j<currentFactorVal.size(); ++j) {
            assert(factors[j](X) == currentFactorVal[j]);
            if(factorToBasisPerturbedValue[j].size() == 2) assert(currentFactorVal[j].second == true);
            infeasibility += !currentFactorVal[j].second;
        }
        assert(infeasibility == currentInfeasibility);
    }


protected:

    // Turns the given constraints into a minimal basis that spans the
    // subspace described by the constraints.
    // Returns a pair whose first element is the basis vectors
    // and whose second element is the dependency mapping from domain element to
    // the set of basis vectors that depend on that element, where the
    // basis vectors are given an index by the ordering of the first element.
    void init(const std::vector<SparseFunction<std::pair<double,bool>,const DOMAIN &>> &distFactors,
              const std::vector<SparseVec<ELEMENT>> &basisVecs, DOMAIN &X) {

        // construct the dependencies from domain index to basis index (i.e. rows of the basis matrix)
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
        for(int factorId = 0; factorId<distFactors.size(); ++factorId) {
            for(int domainIndex: distFactors[factorId].dependencies) {
                for(int basisIndex: domainToBasisDependencies[domainIndex]) {
                    dependencyMatrix.insert(std::pair(factorId, basisIndex));
                }
            }
        }
        std::cout << "mean basis size = " << nBasisEntries * 1.0 / basisVecs.size() << std::endl;
        std::cout << "mean dependency col size = " << dependencyMatrix.size() * 1.0 / basisVecs.size() << std::endl;
        std::cout << "mean dependency row size = " << dependencyMatrix.size() * 1.0/ distFactors.size() << std::endl;
        std::cout << "mean updates per transition = " << pow(dependencyMatrix.size(),2)*1.0/(basisVecs.size()*distFactors.size())  << std::endl;

        // now transfer to forward and backward dependencies
        // set up structure, values will be added later
        factorToBasisPerturbedValue.resize(distFactors.size());
        for(auto &dependencyEntry: dependencyMatrix) {
            factorToBasisPerturbedValue[dependencyEntry.first].emplace_back(dependencyEntry.second);
        }
        basisToFactorPerturbedValue.resize(basisVecs.size());
        for(int factorId = 0; factorId < factorToBasisPerturbedValue.size(); ++factorId) {
            for(RowEntry &entry : factorToBasisPerturbedValue[factorId]) {
                basisToFactorPerturbedValue[entry.basisIndex].emplace_back(factorId, entry.value);
            }
        }

        initState(distFactors,X);
        debug(sanityCheck());
    }

    void initState(const std::vector<SparseFunction<std::pair<double,bool>,const DOMAIN &>> &distFactors, DOMAIN &X) {
        // initialise factor values and infeasibility
        currentFactorVal.reserve(distFactors.size());
        currentInfeasibility = 0;
        for(const auto &factor: distFactors) {
            currentFactorVal.push_back(factor(X)); // Assumes all factors with single basis dependency are feasible
            assert((factorToBasisPerturbedValue[currentFactorVal.size()-1].size() != 2) || (currentFactorVal.back().second == true));
            currentInfeasibility += !currentFactorVal.back().second;
        }

        // initialise weights, ratios and probabilities
        currentWeight.resize(basisVectors.size(), 0.0);
        for(int basisIndex=0; basisIndex < basisVectors.size(); ++basisIndex) {
            X += basisVectors[basisIndex];
            for(ColEntry &colEntry : basisToFactorPerturbedValue[basisIndex]) {
                colEntry.value = distFactors[colEntry.factorIndex](X);
//                std::cout << factorToBasisPerturbedValue[colEntry.factorIndex].size() << std::endl;
                if(factorToBasisPerturbedValue[colEntry.factorIndex].size() == 2 && colEntry.value.second == false) { // if this factor depends only on this basis, stay within feasible region
                    colEntry.value.first -= 24.0; // TODO: test!!!!
                }
                currentWeight[basisIndex] += colEntry.value.first - currentFactorVal[colEntry.factorIndex].first;
            }
            X -= basisVectors[basisIndex];
            basisDistribution.push_back(logWeighttoBasisProb(currentWeight[basisIndex]));
        }

    };

    static double logWeighttoBasisProb(double w) { return w<0.0?exp(w):1.0; } // convert change in energy of a column to probability of proposal


    // Given a transition along basis j,
    // Update each basis, k, for which there exists a factor that depends on both j and k
    void performTransition(int transitionBasis) {
        std::set<int> updatedBasisWeights;
        X += basisVectors[transitionBasis];
        for(ColEntry &colEntry: basisToFactorPerturbedValue[transitionBasis]) {
            int factorIndex = colEntry.factorIndex;
            double newFVal = colEntry.value.first; // factor val after this transition
            double oldFVal = currentFactorVal[factorIndex].first;   // factor val before transition
            for(RowEntry &rowEntry: factorToBasisPerturbedValue[factorIndex]) {
                int basisIndex = rowEntry.basisIndex; // basis of potential perturbation
                std::pair<double,bool> &perturbedFactorVal = rowEntry.value; // factor val after perturbation by basisIndex
                if(basisIndex == transitionBasis) {
                    currentInfeasibility += currentFactorVal[factorIndex].second - perturbedFactorVal.second;
                    currentFactorVal[factorIndex] = perturbedFactorVal;
                    assert((factorToBasisPerturbedValue[factorIndex].size() != 2) || (perturbedFactorVal.second == true));
                }
                currentWeight[basisIndex] += oldFVal - perturbedFactorVal.first; // remove old ratio
                // update the given factor of the given basis transition
                X += basisVectors[basisIndex];
                perturbedFactorVal = factors[factorIndex](X); // record the new perturbed value
                if(factorToBasisPerturbedValue[factorIndex].size() == 2 && perturbedFactorVal.second == false) { // if this factor depends only on this basis (+-), stay within feasible region
                    perturbedFactorVal.first -= 24.0; // TODO: test!!!!
                }
                X -= basisVectors[basisIndex];
                currentWeight[basisIndex] += perturbedFactorVal.first - newFVal; // insert new ratio
                updatedBasisWeights.insert(basisIndex);
            }
        }
        for(int basisIndex: updatedBasisWeights) {
            basisDistribution[basisIndex] = logWeighttoBasisProb(currentWeight[basisIndex]);
        }
    }
};


template<class DISTRIBUTION, class USER,
        class RESULT = std::result_of<USER(std::function<const typename DISTRIBUTION::domain_type &()>)>,
        class SOLVER = decltype(FactorisedDistributionSampler(std::declval<DISTRIBUTION>()))>
static RESULT solve(const DISTRIBUTION &distribution, USER user) {
    FactorisedDistributionSampler sampler(distribution);
    return user(sampler);
}

template<class DISTRIBUTION, class USER,
        class RESULT = std::result_of<USER(std::function<const typename DISTRIBUTION::domain_type &()>)>,
        class SOLVER = decltype(FactorisedDistributionSampler(std::declval<DISTRIBUTION>()))>
static std::future<RESULT> solveAsync(const DISTRIBUTION &distribution, USER user) {
    FactorisedDistributionSampler sampler(distribution);
    return std::async(user,sampler);
}

#endif //ABMCMC_FACTORISEDDISTRIBUTIONSAMPLER_H
