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

template<typename DOMAIN, typename ELEMENT>
class ConstrainedFactorisedSampler {
public:
    typedef std::vector<int> IndexSet;
    typedef std::remove_reference_t<std::remove_const_t<DOMAIN>> stripped_domain;

    std::vector<SparseVec<ELEMENT>> basisVectors;       // the basis vectors by basis index, input index (copy for +-?)
    MutableCategoricalArray         basisDistribution;  // probability of proposing to update the j'th basis
    std::vector<double>             currentLogWeight;   // current log weight of the j'th transition

    std::vector<std::function<std::pair<double,bool>(DOMAIN)>>  factors;
    std::vector<std::pair<double,bool>>                         currentFactorVal;       // current log value of the i'th factor

    std::vector<IndexSet>       dependencyRows; // set of basis vectors that a given factor is dependent on
    std::vector<IndexSet>       dependencyCols; // set of factors that are dependent on a given basis

    stripped_domain             X;      // the current point on the lattice by row.
    int                         currentInfeasibility; // current number of factors that are in widened support
    MCMCStatistics              stats;

    ConstrainedFactorisedSampler<DOMAIN,ELEMENT>(
            ConstrainedFactorisedDistribution<DOMAIN,ELEMENT> &targetDistribution,
            DOMAIN initialSolution):
            X(initialSolution),
            basisDistribution(initialSolution.size() - targetDistribution.constraints.size()),
            dependencyCols(basisDistribution.size()),
            currentInfeasibility(0)
    {
        // initialise factor stuff and dependency matrix
        factors.reserve(targetDistribution.logFactors.size());
        currentFactorVal.reserve(targetDistribution.logFactors.size());
        for(const SparseWidenedFunction<double,DOMAIN> &factor: targetDistribution.logFactors) {
            factors.push_back(factor.widenedFunction);
            currentFactorVal.push_back(factor.widenedValue(X));
            currentInfeasibility += currentFactorVal.back().second;
            dependencyRows.push_back(factor.dependencies);
            for(int dependentBasisIndex: factor.dependencies) {
                dependencyCols[dependentBasisIndex].push_back(factors.size()-1);
            }
        }

        // inisialise basis stuff
        TableauNormMinimiser<ELEMENT> basisMinimiser(targetDistribution);
        basisVectors = basisMinimiser.getMinimalBasisVectors();
        while(basisVectors.size() != basisDistribution.size()) {
            // fill in identity basis for any unconstrained variables at the end of the vector
            basisVectors.push_back(SparseVec<ELEMENT>());
            basisVectors.back().insert(basisVectors.size()-1,1);
        }

        // initialise weights and probabilities
        for(int j=0; j<basisVectors.size(); ++j) {
            currentLogWeight[j] = 0.0;
            X += basisVectors[j];
            for(int i: dependencyCols[j]) {
                currentLogWeight[j] += factors[i](X).first - currentFactorVal[i].first;
            }
            X -= basisVectors[j];
            basisDistribution[j] = logWeighttoBasisProb(currentLogWeight[j]);
        }

    }


    const DOMAIN &nextSample() {
        int proposedBasisIndex = basisDistribution(Random::gen);
        double oldSum = basisDistribution.sum();
        performTransition(proposedBasisIndex);
        double acceptance = oldSum / basisDistribution.sum();
        if(Random::nextDouble() > acceptance) {
            performTransition(proposedBasisIndex ^ 1); // reject: reverse transition
        }
        return X;
    }

protected:
    static double logWeighttoBasisProb(double w) { return w<0.0?exp(w):1.0; } // convert change in energy of a column to probability of proposal


    // Given a transition along basis j,
    // Update each basis, k, for which there exists a factor that depends on both j and k
    void performTransition(int perturbedBasis) {
        X += basisVectors[perturbedBasis];
        for(int factorIndex: dependencyCols[perturbedBasis]) {
            for(int basisIndex: dependencyRows[factorIndex]) {
                // update the given factor of the given basis transition
                std::pair<double,bool> oldLogF = currentFactorVal[factorIndex];
                currentFactorVal[factorIndex] = factors[factorIndex](X);
                X += basisVectors[basisIndex];
                double perturbedLogFVal = factors[factorIndex](X).first;
                X -= basisVectors[basisIndex];
                currentLogWeight[basisIndex] += oldLogF.first - 2.0 * currentFactorVal[factorIndex].first + perturbedLogFVal;
                basisDistribution[basisIndex] = logWeighttoBasisProb(currentLogWeight[basisIndex]);
                currentInfeasibility += currentFactorVal[factorIndex].second - oldLogF.second;
            }
        }
    }
};


#endif //ABMCMC_CONSTRAINEDFACTORISEDSAMPLER_H
