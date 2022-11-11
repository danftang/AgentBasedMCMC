// Samples from an ABM posterior by using an exponential proposal function in the form
// P'(B) = A\prod_i e^{w_iB_i}
// subject to
// MB >= 0 (i.e. X_i >= 0 for all i)
// where
// B is the coordinate in basis-vector space and M is a matrix whose columns are the basis vectors.
// So that the probability of proposing a transition from B to B+1_i (i.e. incrementing the i^th basis vector)
//
// P(B->B+1_i) = \sqrt{P'(B+1_i)/P'(B)} / s(B)
// where
// s(B) = \sum_j \sqrt{P'(B+1_j)/P'(B)}
//
// In order to ensure mixing, the sampling process is a hybrid of Metropolis-Hastings and rejection sampling.
// We perform MH on max(P'(B),P(B)) then rejection where P'(B) > P(B). So, acceptance
//
// \alpha(B->B+1_i) = max(P'(B+1_i),P(B+1_i))P(B+1_i -> B) / max(P'(B),P(B))P(B -> B+1_i)
//
//   = (P'(B) / max(P'(B),P(B))) (max(P'(B+1_i),P(B+1_i)) / P'(B+1_i)) (s(B+1_j) / s(B))
//
// If accepted, the Markov state is updated, however, a token is only output with probability
// r = P(B) / max(P'(B),P(B))
//
// Created by daniel on 11/11/22.
//

#ifndef ABMCMC_ABMEXPSAMPLER_H
#define ABMCMC_ABMEXPSAMPLER_H

#include <vector>
#include <functional>

#include "include/SparseVec.h"
#include "include/MutableCategoricalArray.h"
#include "diagnostics/MultiChainStats.h"

template<typename DOMAIN, typename ELEMENT = typename DOMAIN::value_type>
class ABMExpSampler {
public:

    // USER is a function object that takes a producer of samples and returns a result
//    template<class USER, class RESULT = typename SamplingProblem<ConstrainedFactorisedDistribution<DOMAIN,ELEMENT>,USER>::result_type>
//    static RESULT solve(const SamplingProblem<ConstrainedFactorisedDistribution<DOMAIN,ELEMENT>,USER> &problem) {
//        FactorisedDistributionSampler sampler(problem.distribution);
//        return problem.user(sampler);
//    }



public:
    typedef std::vector<int> IndexSet;

    std::vector<SparseVec<ELEMENT>> basisVectors;       // the basis vectors by basis index and domain index (doubled for +-)
    std::vector<std::function<double(const DOMAIN &)>>  factors; // Factorised target distribution
    std::vector<double>  currentFactorVal;  // current log value of the i'th factor
    std::vector<IndexSet>     basisToFactorDependencies; // set of factors that are dependent on a given basis

    MutableCategoricalArray         basisDistribution;  // probability of proposing to update the j'th basis
    std::vector<double>             basisWeights;   // exp coefficients of the proposal distribution
    std::vector<int>                basisInfeasibility; // number of -ve entries in a basis that would currently cause infeasibility
    std::vector<IndexSet>           domainToNegativeBasisEntry; // basis vectors that have a negative entry for a given domain index.

    DOMAIN                      X;      // the current Markov state
    MCMCStatistics              stats;  // statistics on all samples since construction

public:
    ABMExpSampler(ABMExpSampler &&other)=delete;
    ABMExpSampler(const ABMExpSampler &other)=delete;

    template<class STARTSTATE>
    ABMExpSampler(const ABMPosterior<DOMAIN,STARTSTATE> &targetDistribution):
            ABMExpSampler(targetDistribution.basis, targetDistribution.factors) { }

    template<class STARTSTATE>
    ABMExpSampler(const ABMPrior<DOMAIN,STARTSTATE> &targetPrior, double constraintKappa = 6):
            ABMExpSampler(Basis(targetPrior.approximateEntropiesByVarIndex(),targetPrior, constraintKappa), targetPrior.factors) { }

    ABMExpSampler(
            const Basis<DOMAIN> &basis,
            const std::vector<double> &basisWeights,
            const std::vector<SparseFunction<double,const DOMAIN &>> &distFactors): basisWeights(basisWeights), X(basis.origin) {

        // Choose random start state on the basis manifold (X is already at basis origin)
        Random::gen.seed(Random::nextRandomSeed());
        for(int i=0; i<basis.basisVectors.size(); ++i) {
            if(Random::nextBool(0.005)) X += basis.basisVectors[i];
        }

        // double-up unique basis vectors for + and -
        basisVectors.reserve(basis.basisVectors.size()*2);
        for(int i=0; i < basis.basisVectors.size(); ++i) {
            basisVectors.push_back(-(basis.basisVectors[i]));
            basisVectors.push_back( basis.basisVectors[i]);
        }

        // copy factors
        factors.reserve(distFactors.size());
        for(const auto &factor: distFactors) factors.push_back(factor);

        // construct the dependencies from domain index to basis index (i.e. rows of the basis matrix)
        std::vector<IndexSet> domainToBasisDependencies(DOMAIN::dimension);
        int nBasisEntries = 0;

        for(int basisIndex=0; basisIndex < basisVectors.size(); ++basisIndex) {
            for(const auto &[domainIndex, entryVal]: basisVectors[basisIndex]) {
                domainToBasisDependencies[domainIndex].push_back(basisIndex);
                if(entryVal < 0) domainToNegativeBasisEntry[domainIndex].push_back(basisIndex);
                ++nBasisEntries;
            }
        }

//        // calculate factor by basis-dependency matrix
//        std::set<std::pair<int,int>> dependencyMatrix; // entries
//        for(int factorId = 0; factorId<distFactors.size(); ++factorId) {
//            for(int domainIndex: distFactors[factorId].dependencies) {
//                for(int basisIndex: domainToBasisDependencies[domainIndex]) {
//                    dependencyMatrix.insert(std::pair(factorId, basisIndex));
//                }
//            }
//        }

        // initialise basisToFactorDependencies
        basisToFactorDependencies.resize(basisVectors.size());
        for(int factorId = 0; factorId < factors.size(); ++factorId) {
            for(int domainIndex : distFactors[factorId].dependencies) {
                for(int basisIndex: domainToBasisDependencies[domainIndex]) {
                    basisToFactorDependencies[basisIndex].emplace_back(factorId);
                }
            }
        }
        for(int basisIndex = 0; basisIndex<basisVectors.size(); ++basisIndex) {
            std::sort(basisToFactorDependencies[basisIndex].begin(), basisToFactorDependencies[basisIndex].end());
            auto newEnd = std::unique(basisToFactorDependencies[basisIndex].begin(), basisToFactorDependencies[basisIndex].end());
            basisToFactorDependencies[basisIndex].resize(newEnd - basisToFactorDependencies[basisIndex].begin());
        }

        // initialise factor values
        currentFactorVal.reserve(distFactors.size());
        for(const auto &factor: distFactors) currentFactorVal.push_back(factor(X));

        // initialise basis infeasibilities
        basisInfeasibility.resize(basisVectors.size(), 0);
        for(int domainIndex=0; domainIndex < DOMAIN::size(); ++domainIndex) {
            if(X[domainIndex] == 0) {
                for(int basisIndex: domainToNegativeBasisEntry[domainIndex]) {
                    ++basisInfeasibility[basisIndex];
                }
            }
        }

        // initialise basis probability distribution
        for(int basisIndex=0; basisIndex < basisVectors.size(); ++basisIndex) {
            double basisProb = (basisInfeasibility[basisIndex]==0?exp(basisWeights[basisIndex]):0.0);
            basisDistribution.push_back(basisProb);
        }
        debug(sanityCheck());

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

    double currentLogProb() {
        double logProb = 0.0;
        for(int i=0; i<factors.size(); ++i) {
            auto factorVal = factors[i](X);
            logProb += factorVal.first;
            if(factorVal.second == false && factorToBasisPerturbedValue[i].size() == 2 && currentFactorVal[i].second == true) {
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
            std::pair<double,bool> realFactorVal = factors[j](X);
            if(realFactorVal.second == false && factorToBasisPerturbedValue[j].size() == 2 && currentFactorVal[j].second == true) {
                assert(currentFactorVal[j].first == realFactorVal.first - 24.0);
                assert(currentFactorVal[j].second == false);
            } else {
                assert(currentFactorVal[j] == realFactorVal);
            }
            infeasibility += !realFactorVal.second;
        }
        assert(infeasibility == currentInfeasibility);
    }


protected:


    // convert the log(P(destination)/P(source)) of a transition to the probability of proposal
    static double logWeighttoBasisProb(double w) {
        return w<0.0?exp(w):1.0;
//        return exp(0.5*w); // square root of prob ratio
    }

    bool basisIsActive(int basisIndex) {

    }


    // Given a transition along basis j,
    // Update each basis, k, for which there exists a factor that depends on both j and k
    void performTransition(int transitionBasis) {
        std::set<int> updatedBasisWeights;
        X += basisVectors[transitionBasis];
        for(ColEntry &colEntry: basisToFactorDependencies[transitionBasis]) {
            int factorIndex = colEntry.factorIndex;
            double newFVal = colEntry.value.first; // factor val after this transition
            double oldFVal = currentFactorVal[factorIndex].first;   // factor val before transition
            for(RowEntry &rowEntry: factorToBasisPerturbedValue[factorIndex]) {
                int basisIndex = rowEntry.basisIndex; // basis of potential perturbation
                std::pair<double,bool> &perturbedFactorVal = rowEntry.value; // factor val after perturbation by basisIndex
                if(basisIndex == transitionBasis) {
                    currentInfeasibility += currentFactorVal[factorIndex].second - perturbedFactorVal.second;
                    currentFactorVal[factorIndex] = perturbedFactorVal;
                    assert(!(perturbedFactorVal.second == false && factorToBasisPerturbedValue[factorIndex].size() == 2 && currentFactorVal[factorIndex].second == true));
                }
                currentWeight[basisIndex] += oldFVal - perturbedFactorVal.first; // remove old ratio
                // update the given factor of the given basis transition
                X += basisVectors[basisIndex];
                perturbedFactorVal = factors[factorIndex](X); // record the new perturbed value
                if(perturbedFactorVal.second == false && factorToBasisPerturbedValue[factorIndex].size() == 2 && currentFactorVal[factorIndex].second == true) {
                    // if this factor depends only on this basis (i.e. + and - basis vectors), don't transition from feasible to infeasible
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



#endif //ABMCMC_ABMEXPSAMPLER_H
