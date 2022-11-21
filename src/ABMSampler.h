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

#ifndef ABMCMC_ABMSAMPLER_H
#define ABMCMC_ABMSAMPLER_H

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
class ABMSampler {
public:

    // USER is a function object that takes a producer of samples and returns a result
//    template<class USER, class RESULT = typename SamplingProblem<ConstrainedFactorisedDistribution<DOMAIN,ELEMENT>,USER>::result_type>
//    static RESULT solve(const SamplingProblem<ConstrainedFactorisedDistribution<DOMAIN,ELEMENT>,USER> &problem) {
//        FactorisedDistributionSampler sampler(problem.distribution);
//        return problem.user(sampler);
//    }



public:
    typedef std::vector<int> IndexSet;

    // Rows and cols of the factor/basis dependency matrix
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
    ABMSampler(ABMSampler &&other)=delete;
    ABMSampler(const ABMSampler &other)=delete;
//    basisVectors(other.basisVectors),
//    factors(other.factors),
//    X(other.X) {
//        init(factors, basisVectors, X);
//    }

//    ABMPosteriorSampler(const ConstrainedFactorisedDistribution<DOMAIN, ELEMENT> &targetDistribution):
//            ABMPosteriorSampler(Basis(targetDistribution), targetDistribution.factors) { }
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
    ABMSampler(const ABMPosterior<DOMAIN,STARTSTATE> &targetDistribution):
            ABMSampler(targetDistribution.basis, targetDistribution.factors) { }

    template<class STARTSTATE>
    ABMSampler(const ABMPrior<DOMAIN,STARTSTATE> &targetPrior, double constraintKappa = 6):
            ABMSampler(Basis(targetPrior.approximateEntropiesByVarIndex(),targetPrior, constraintKappa), targetPrior.factors) { }


    ABMSampler(
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
        std::cout << "Finding initial feasible solution..." << std::endl;
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
//            std::cout << "Taking sample. Current infeasibility = " << currentInfeasibility << std::endl;
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
//            debug(sanityCheck());
        } while(currentInfeasibility != 0);
        return X;
    }


    // Log prob of current state X
    double currentLogProb() {
        double logProb = 0.0;
        for(int i=0; i<factors.size(); ++i) {
            auto factorVal = factors[i](X);
            logProb += factorVal.first;
//            if(factorVal.second == false && factorToBasisPerturbedValue[i].size() == 2 && currentFactorVal[i].second == true) {
//                logProb -= 24.0; // TODO: Test!!!!
//            }
        }
        return logProb;
    }



    void sanityCheck() {
        // check weights
        double logProbX0 = currentLogProb();
        for(int basisId=0; basisId < basisVectors.size(); ++basisId) {
            X += basisVectors[basisId];
            double logRatioOfProbs = currentLogProb() - logProbX0;
            assert(fabs(currentWeight[basisId] - logRatioOfProbs) < 1e-8);
            X -= basisVectors[basisId];
            assert(fabs(basisProbabilityFromCurrentWeight(basisId) - basisDistribution[basisId]) < 1e-8);
        }
        int infeasibility = 0;
        for(int j=0; j<currentFactorVal.size(); ++j) {
            std::pair<double,bool> realFactorVal = factors[j](X);
//            if(realFactorVal.second == false && factorToBasisPerturbedValue[j].size() == 2 && currentFactorVal[j].second == true) {
//                assert(currentFactorVal[j].first == realFactorVal.first - 24.0);
//                assert(currentFactorVal[j].second == false);
//            } else {
                assert(currentFactorVal[j] == realFactorVal);
//            }
            infeasibility += !realFactorVal.second;
        }
        assert(infeasibility == currentInfeasibility);
    }


protected:

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
            currentFactorVal.push_back(factor(X));
            // Assumes all factors with single basis dependency are feasible
//            assert((factorToBasisPerturbedValue[currentFactorVal.size()-1].size() != 2) || (currentFactorVal.back().second == true));
            currentInfeasibility += !currentFactorVal.back().second;
        }

        // initialise weights, ratios and probabilities
        currentWeight.resize(basisVectors.size(), 0.0);
        for(int basisIndex=0; basisIndex < basisVectors.size(); ++basisIndex) {
            X += basisVectors[basisIndex];
            for(ColEntry &colEntry : basisToFactorPerturbedValue[basisIndex]) {
                colEntry.value = distFactors[colEntry.factorIndex](X);
//                std::cout << factorToBasisPerturbedValue[colEntry.factorIndex].size() << std::endl;
//                if(colEntry.value.second == false && factorToBasisPerturbedValue[colEntry.factorIndex].size() == 2 && currentFactorVal[colEntry.factorIndex].second == true) {
//                    // if this factor depends only on this basis (i.e. + and - basis vectors), don't transition from feasible to infeasible
//                    colEntry.value.first -= 24.0; // TODO: test!!!!
//                }
                currentWeight[basisIndex] += colEntry.value.first - currentFactorVal[colEntry.factorIndex].first;
            }
            X -= basisVectors[basisIndex];
            basisDistribution.push_back(basisProbabilityFromCurrentWeight(basisIndex));
        }

    };


    double basisProbabilityFromCurrentWeight(int basisIndex) {
        ELEMENT XiNonBasic = X[basisVectors[basisIndex].indices[0]];
        return ((XiNonBasic < 0) || (XiNonBasic + basisVectors[basisIndex].values[0] >= 0))?
            logWeighttoBasisProb(currentWeight[basisIndex]):0.0;
    }

    // convert the log(P(destination)/P(source)) of a transition to the probability of proposal
    static double logWeighttoBasisProb(double w) {
//        return w<0.0?exp(w):1.0;
        return exp(0.5*w); // square root of prob ratio
    }


    // Given a transition along basis j,
    // Update each basis, k, for which there exists a factor that depends on both j and k
    void performTransition(int transitionBasis) {
        std::vector<int> updatedBasisIndices;
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
//                    assert(!(perturbedFactorVal.second == false && factorToBasisPerturbedValue[factorIndex].size() == 2 && currentFactorVal[factorIndex].second == true));
                }
                currentWeight[basisIndex] += oldFVal - perturbedFactorVal.first; // remove old ratio
                // update the given factor of the given basis transition
                X += basisVectors[basisIndex];
                perturbedFactorVal = factors[factorIndex](X); // record the new perturbed value
//                if(perturbedFactorVal.second == false && factorToBasisPerturbedValue[factorIndex].size() == 2 && currentFactorVal[factorIndex].second == true) {
//                    // if this factor depends only on this basis (i.e. + and - basis vectors), don't transition from feasible to infeasible
//                    perturbedFactorVal.first -= 24.0; // TODO: test!!!!
//                }
                X -= basisVectors[basisIndex];
                currentWeight[basisIndex] += perturbedFactorVal.first - newFVal; // insert new ratio
                updatedBasisIndices.push_back(basisIndex);
            }
        }
        std::sort(updatedBasisIndices.begin(), updatedBasisIndices.end());
        auto newEnd = std::unique(updatedBasisIndices.begin(), updatedBasisIndices.end());
        for(auto it=updatedBasisIndices.begin(); it != newEnd; ++it) {
            int basisIndex = *it;
            basisDistribution[basisIndex] = basisProbabilityFromCurrentWeight(basisIndex);
        }
    }
};


template<class DISTRIBUTION, class USER,
        class RESULT = std::result_of<USER(std::function<const typename DISTRIBUTION::domain_type &()>)>,
        class SOLVER = decltype(ABMSampler(std::declval<DISTRIBUTION>()))>
static RESULT solve(const DISTRIBUTION &distribution, USER user) {
    ABMSampler sampler(distribution);
    return user(sampler);
}

template<class DISTRIBUTION, class USER,
        class RESULT = std::result_of<USER(std::function<const typename DISTRIBUTION::domain_type &()>)>,
        class SOLVER = decltype(ABMSampler(std::declval<DISTRIBUTION>()))>
static std::future<RESULT> solveAsync(const DISTRIBUTION &distribution, USER user) {
    ABMSampler sampler(distribution);
    return std::async(user,sampler);
}

#endif //ABMCMC_ABMSAMPLER_H
