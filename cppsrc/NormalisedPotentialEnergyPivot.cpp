//
// Created by daniel on 06/07/2021.
//

#include <cfloat>

#include "Random.h"
#include "StlStream.h"
#include "NormalisedPotentialEnergyPivot.h"
#include "SimplexMCMC.h"
#include "agents/PredPreyAgent.h"

NormalisedPotentialEnergyPivot::NormalisedPotentialEnergyPivot(SimplexMCMC &simplex)
:
ProposalPivot(simplex,0,0),
//kappaCol(std::max(std::log(p1*simplex.nNonBasic()),1.0)),  // kappaCol(3.5) {
//kappaCol(std::max(std::log((simplex.nNonBasic()-1.0)/(1.0/p1 - 1.0)),1.0)),  // kappaCol(3.5) {
cdf(simplex.nNonBasic() + 1)
{
    // setup initial cached values
    initCache();
}

void NormalisedPotentialEnergyPivot::initCache() {
    col.resize(simplex.nBasic()+1);
    deltaj = 0.0;
    infeasibilityCosts.emplace_front(simplex.nBasic()+1);
    calcProposedInfeasibilityCosts();
    if(infeasibilityCosts.size() > 1) infeasibilityCosts.pop_front();
    calcProposedReducedCosts();
    calcProposedEnergies();
    cdf[0] = 0.0;
    recalculateCDF();
}

const NormalisedPotentialEnergyPivot &NormalisedPotentialEnergyPivot::nextProposal() {
    // update cache
    if(simplex.lastSampleWasAccepted) {
        // remove last state's cache
        if(infeasibilityCosts.size() > 1)   infeasibilityCosts.pop_front();
        if(reducedCosts.size() > 1) reducedCosts.pop_front();
        if(potentialEnergies.size() > 1) {
            potentialEnergies.pop_front();
            recalculateCDF();
        }
//        if(totalPotentials.size() >1) totalPotentials.pop_front();
        if(energyNormalisation.size() >1) energyNormalisation.pop_front();
    } else {
        // remove cache for last proposal
        if(infeasibilityCosts.size() > 1)   infeasibilityCosts.pop_back();
        if(reducedCosts.size() > 1) reducedCosts.pop_back();
        if(potentialEnergies.size() > 1)    potentialEnergies.pop_back();
//        if(totalPotentials.size() >1) totalPotentials.pop_back();
        if(energyNormalisation.size() >1) energyNormalisation.pop_back();
    }
//    checkCurrentCacheIsValid();
    chooseCol();
    chooseRow();
    calcAcceptanceContrib();
    return *this;
}


void NormalisedPotentialEnergyPivot::recalculateCDF() {
    double cumulativeP = 0.0;
    const std::vector<double> &currentPotentials = potentialEnergies.front();
//    int basisEnergy = 0;
    for(int q=1; q <= simplex.nNonBasic(); ++q) {
//        double colPotential = potentialEnergy(q, currentReducedCost);
//        Ep += colPotential;
//        int isBasis = isInPredPreyPreferredBasis(q);
//        if(isBasis) ++basisEnergy;
        cumulativeP += exp(kappaCol* currentPotentials[q]);
        cdf[q] = cumulativeP;
    }
//    std::cout << "Basis energy = " << basisEnergy << std::endl;
//    assert(cumulativeP > simplex.nNonBasic()); // there should be at least one high energy column
}

bool NormalisedPotentialEnergyPivot::isInPredPreyPreferredBasis(int j) {
    int kProb = simplex.kSimTokProb[j+simplex.nNonBasic()];
    return kProb>simplex.nBasic() && kProb%PredPreyAgentBase::actDomainSize() == PredPreyAgentBase::DIE;
}


// Chooses columns based on the exponential of their "Potential energy" which is defined as sign(d_j)*(2X_j-1)
// where d_j is the reduced cost of the column and X_j is its value, and a column is chosen with
// a probaqbility equal to e^kEj/sum_l e^kEl, where Ej is the potential energy of column j.
//
// The porobability of a state is penalised by A * sum_l e^-kEl * e^-kE / n where E is the total potential energy
// summed over all non-basic columns and n is the number of non-basics. E is observed to be non-negative
// (TODO: should we limit penalty to be always less than 1?).
//
// As a consequence, the contribution of the column choice to the acceptance probability is
// e^-k(DE - DEj) where DE and DEj are the change in total and column energy respectively, during the pivot.
void NormalisedPotentialEnergyPivot::chooseCol() {
    auto chosenJ = std::lower_bound(
            cdf.begin(),
            cdf.end(),
            Random::nextDouble(0.0, cdf[simplex.nNonBasic()])
        );
    setCol(chosenJ - cdf.begin());
}

void NormalisedPotentialEnergyPivot::chooseRow() {
    std::multimap<double, int> transitions = getPivotsByDeltaJ(); // maps from delta_j to LogPMF-index.

//    std::cout << "Transition range = " << transitions.begin()->first << " to " << transitions.rbegin()->first << std::endl;

    int lowestDeltaJ = transitions.begin()->first;
//    int highestDeltaJ = transitions.rbegin()->first;

    // now populate pivotPMF by going from lowest to highest delta_j in order
    std::vector<double> pivotPMF(nonZeroRows.size() * 2 + 2,
                                 DBL_MAX); // index is (2*nonZeroRowIndex + toUpperBound), final two are lower and upper bound swap, value is probability mass
    std::map<double, double> reducedCostByDj;//( highestDeltaJ - lowestDeltaJ + 1, 0.0);
    double lastDj = lowestDeltaJ;
    double infeas = 0.0;
    double infeasMin = DBL_MAX;
//    double dDf_dDj = colInfeasibilityGradient(lastDj - 2.0*tol);
    double DdDf_dDj;
    double upperRateOfChange = 0.0; // rate of change stored in transitions from feasible to infeasible (as opposed to infeasible to feasible)
    reducedCostByDj[lowestDeltaJ] = colInfeasibilityGradient(lastDj - 2.0 * tol);
    for (auto[Dj, pmfIndex]: transitions) {
        double pivotPointVal = 1.0;
        if (pmfIndex < nonZeroRows.size() * 2) {
            pivotPointVal = col[nonZeroRows[pmfIndex / 2]];
            DdDf_dDj = fabs(pivotPointVal);
        } else {
            // null pivot or bound swap
            DdDf_dDj = 1.0;
        }
        if (Dj != lastDj) {
            reducedCostByDj[Dj] = reducedCostByDj[lastDj] + upperRateOfChange;
            infeas += (Dj - lastDj) * reducedCostByDj[Dj];
            lastDj = Dj;
            upperRateOfChange = 0.0;
        }
        bool varAtUpperBound = pmfIndex % 2;
        if ((pivotPointVal < 0.0) xor varAtUpperBound) {
            // is an upper var
            upperRateOfChange += DdDf_dDj;
        } else {
            //is a lower var
            reducedCostByDj[Dj] += DdDf_dDj;
        }

        if (isActive(pmfIndex)) {
            pivotPMF[pmfIndex] = infeas;
            if (infeas < infeasMin) infeasMin = infeas;
        }
//        std::cout << "deltaj = " << Dj << " pivotProb = " << pivotPMF[pmfIndex] << " infeas = " << infeas << " dDf_dDj = " << reducedCostByDj[Dj] << std::endl;
    }
    assert(infeasMin != DBL_MAX);

    // now add proposed j^th col potential to active row probabilities
    for (auto[Dj, pmfIndex]: transitions) {
        if (pivotPMF[pmfIndex] < DBL_MAX) {
            double pivotPointVal = (pmfIndex < nonZeroRows.size() * 2) ? col[nonZeroRows[pmfIndex / 2]] : 1.0;
            double proposedReducedCost = reducedCostByDj[Dj] * pivotPointVal;
            double proposedPotential = potentialEnergy(pmfIndex % 2, proposedReducedCost);
            pivotPMF[pmfIndex] += (kappaCol / kappaRow) * proposedPotential;
        }
    }

    // rescale, take exponential and turn pivotPMF into a cumulative distribution
    double cumulativeMass = 0.0;
    for(int i=0; i<pivotPMF.size(); ++i) {
        cumulativeMass += exp(kappaRow*(pivotPMF[i] - infeasMin));
        pivotPMF[i] = cumulativeMass;
    }

    int chosenRowIndex = std::lower_bound(
            pivotPMF.begin(),
            pivotPMF.end(),
            Random::nextDouble(0.0, pivotPMF.back())) - pivotPMF.begin();
    setToPivotIndex(chosenRowIndex);

}

//void NormalisedPotentialEnergyPivot::chooseRow() {
//    std::multimap<double, int> transitions = getPivotsByDeltaJ(); // maps from delta_j to LogPMF-index.
//
//    // now populate pivotPMF by going from lowest to highest delta_j in order
//    std::vector<double> pivotPMF(nonZeroRows.size() * 2 + 2, DBL_MAX); // index is (2*nonZeroRowIndex + toUpperBound), final two are lower and upper bound swap, value is probability mass
//    double lastDj = transitions.begin()->first;
//    double infeas = 0.0; //infeasibility(lastDj) - infeasibility(0.0);
//    double infeasMin = DBL_MAX;
//    double dDf_dDj = colInfeasibilityGradient(lastDj - 2.0*tol);
//    for(auto [Dj, pmfIndex] : transitions) {
//        infeas += (Dj - lastDj) * dDf_dDj;
//        lastDj = Dj;
//        if(isActive(pmfIndex)) {
//            pivotPMF[pmfIndex] = infeas;
//            if(infeas < infeasMin) infeasMin = infeas;
//        }
////        std::cout << "deltaj = " << Dj << " pivotProb = " << pivotPMF[pmfIndex] << " infeas = " << infeas << " dDf_dDj = " << dDf_dDj << " Feasibility = " << infeasibility(Dj) << std::endl;
//        if(pmfIndex < nonZeroRows.size()*2) {
//            dDf_dDj += fabs(col[nonZeroRows[pmfIndex / 2]]);
//        } else {
//            dDf_dDj += 1.0;
//        }
//    }
//    assert(infeasMin != DBL_MAX);
//
//    // rescale and take exponential
//    for(int i=0; i<pivotPMF.size(); ++i) {
//        pivotPMF[i] = exp(kappaRow*(pivotPMF[i] - infeasMin));
//    }
////    pivotPMF[2*nonZeroRows.size() + simplex.isAtUpperBound(j)] *= 0.01; // preference against null pivot if other possibilities exist
//
//    // nextIntFromDiscrete row
//
//    setToPivotIndex(Random::nextIntFromDiscrete(pivotPMF.begin(), pivotPMF.end()));
//
//    //    std::cout << "Chose pivot with prob " << pivotPMF[pivotChoice] << " Delta = " << deltaj << std::endl;
////  std::cout << "Proposing pivot from " << !sourceObjectiveIsZero << " to " << destinationObjective << std::endl;
//}



// Calculates log acceptance contribution which is kappaCol times
// the proposed change in the potential energy, ignoring column j.
void NormalisedPotentialEnergyPivot::calcAcceptanceContrib() {

    // if deltaj != 0.0, then the infeasibility of the variables may change so recalculate
    // infeasibilityCost of proposal
    if(deltaj != 0.0) calcProposedInfeasibilityCosts();

    // need to recalculate potential energy if
    //   - infeasibility cost has changed
    //   - we're doing a real pivot on a column with non-zero reduced cost
    if(infeasibilityCosts.size() > 1 || (i>0 && reducedCosts.front()[j] != 0.0)) {
        calcProposedReducedCosts();
        calcProposedEnergies();
    }

    // logAcceptanceContribution is defined as -k times the change in total potential energy in all cols q!=j
    if(potentialEnergies.size() > 1) {
//        logAcceptanceContribution = kappaCol * (
//                (totalPotentials.front() - potentialEnergies.front()[j]) -
//                (totalPotentials.back()  - potentialEnergies.back()[j])
//                );

        logAcceptanceContribution = log(energyNormalisation.front()/energyNormalisation.back());

//        std::cout << "normalisation ratio = " << energyNormalisation.front()/energyNormalisation.back() << std::endl;

    } else {
        logAcceptanceContribution = 0.0;
    }

}


// Calculates the post-pivot infeasibility costs for each row and pushes to
// the cache.
void NormalisedPotentialEnergyPivot::calcProposedInfeasibilityCosts() {
    bool hasChanged = false;
    const std::vector<double> &currentInfeasibilityCost = infeasibilityCosts.front();
    std::vector<double> proposedInfeasibilityCost(simplex.nBasic()+1);
    proposedInfeasibilityCost[0] = 0.0;
    for (int p = 1; p <= simplex.nBasic(); ++p) {
        int k;
        double postPivotBetap;
        if (p != i) {
            k = simplex.head[p];
            postPivotBetap = simplex.beta[p] + deltaj * col[p];
        } else {
            k = simplex.head[simplex.nBasic() + j];
            postPivotBetap = (simplex.isAtUpperBound(j) ? simplex.u[k] : simplex.l[k]) + deltaj;
        }
//        proposedInfeasibilityCost[p] = infeasibilityGradient(postPivotBetap, simplex.l[k], simplex.u[k]);
        proposedInfeasibilityCost[p] = infeasibilityCostFn(postPivotBetap, simplex.l[k], simplex.u[k]);
        if(currentInfeasibilityCost[p] != proposedInfeasibilityCost[p]) hasChanged = true;
    }
    if(hasChanged) infeasibilityCosts.push_back(std::move(proposedInfeasibilityCost));
}


// Post-pivot reduced costs are defined as
// R' = -C'B'^-1N,
// where
// C' is the post-pivot infeasibility cost,
// B'^-1 is the inverse of the post-pivot basis matrix
// and N is the non-basic coefficient matrix
// We can express C'B'^-1 as (C'E)B^-1 where B^-1 is the pre-pivot basis inverse
// and E is the elementary transform matrix EB^-1 = B'^-1, which is the identity
// matrix in all but the i'th column.
void NormalisedPotentialEnergyPivot::calcProposedReducedCosts() {
    std::vector<double> proposedCBI(infeasibilityCosts.back()); // will eventually be proposed C'B'^-1
    if (i > 0) {
        // if we're actually pivoting,
        // proposedCBI = proposedCBI * E, where E is the elementary transform matrix (see above).
        double CEi = 0.0;
        for (int p = 1; p <= simplex.nBasic(); ++p) {
            if (p != i) {
                CEi -= proposedCBI[p] * col[p] / col[i];
            } else {
                CEi -= proposedCBI[p] / col[i]; // -ve because col is tableau, not B^-1N_j
            }
        }
        proposedCBI[i] = CEi;

        simplex.btran(proposedCBI); // in-place post-multiplication by B^-1
        std::swap(simplex.head[i], simplex.head[simplex.m+j]);
        std::vector<double> proposedReducedCost = simplex.piTimesMinusN(proposedCBI);
        std::swap(simplex.head[i], simplex.head[simplex.m+j]);
        // now add feasible energy cost of each non-basic
//        for(int j=1; j<proposedReducedCost.size(); ++j) proposedReducedCost[j] += feasibleEnergy;
        reducedCosts.push_back(std::move(proposedReducedCost));
    } else {
        simplex.btran(proposedCBI); // in-place post-multiplication by B^-1
        std::vector<double> proposedReducedCost = simplex.piTimesMinusN(proposedCBI);
        reducedCosts.push_back(std::move(proposedReducedCost));
    }
//    std::cout << "Proposed reduced costs : " << reducedCosts.back() << std::endl;
}


// Given that proposed reduced costs have been calculated, calculate the
// proposed potential energies and total energy and push to the cache
void NormalisedPotentialEnergyPivot::calcProposedEnergies() {
//    double proposedEp = 0.0;
    double proposedEnergyNorm = 0.0;
    std::vector<double> proposedPotential(simplex.nNonBasic()+1);
    std::vector<double> &proposedReducedCost = reducedCosts.back();
    for (int q = 1; q <= simplex.nNonBasic(); ++q) {
        bool proposedQIsAtUpperBound = (q == j)?leavingVarToUpperBound:simplex.isAtUpperBound(q);
        double E = proposedPotential[q] = potentialEnergy(proposedQIsAtUpperBound, proposedReducedCost[q]);
//        proposedEp += E;
        proposedEnergyNorm += exp(kappaCol*E);
    }
    potentialEnergies.push_back(std::move(proposedPotential));
//    totalPotentials.push_back(proposedEp);
    energyNormalisation.push_back(proposedEnergyNorm);
}

// Sanity check for calcAcceptanceContrib()
void NormalisedPotentialEnergyPivot::calcAcceptanceContribCheck() {

    bool jIsAtUpperBound = simplex.isAtUpperBound(j);
    simplex.pivot(*this); // pivot(i, j, col, leavingVarToUpperBound);
    infeasibilityCosts.push_back(infeasibilityCost());
    std::vector<double> proposedCost = infeasibilityCosts.back();
    simplex.btran(proposedCost);
    std::vector<double> proposedReducedCost = simplex.piTimesMinusN(proposedCost);
    double postPivotEp = 0.0;
    std::vector<double> proposedPotentialEnergy(simplex.nNonBasic()+1);
    proposedPotentialEnergy[0] = 0.0;
    for(int q=1; q <= simplex.nNonBasic(); ++q) {
        proposedPotentialEnergy[q] = potentialEnergy(simplex.isAtUpperBound(q),proposedReducedCost[q]);
        postPivotEp += proposedPotentialEnergy[q];
    }
    reducedCosts.push_back(std::move(proposedReducedCost));
    potentialEnergies.push_back(std::move(proposedPotentialEnergy));
//    totalPotentials.push_back(postPivotEp);
    std::vector<double> revCol = simplex.tableauCol(j);
    simplex.pivot(i, j, revCol, jIsAtUpperBound);

    // TODO: update this for energyNormalisation
//    logAcceptanceContribution = kappaCol*(
//            (totalPotentials.front() - potentialEnergies.front()[j]) -
//            (totalPotentials.back() - potentialEnergies.back()[j]));
}

void NormalisedPotentialEnergyPivot::checkCurrentCacheIsValid() {
    std::vector<double> currentCost = infeasibilityCost();
    assert(currentCost == infeasibilityCosts.front());
    simplex.btran(currentCost);
    std::vector<double> currentReducedCost = simplex.piTimesMinusN(currentCost);
    assert(currentReducedCost == reducedCosts.front());
    double Enorm = 0.0;
    std::vector<double> currentPotentialEnergy(simplex.nNonBasic()+1);
    currentPotentialEnergy[0] = 0.0;
    for(int q=1; q <= simplex.nNonBasic(); ++q) {
        currentPotentialEnergy[q] = potentialEnergy(simplex.isAtUpperBound(q),currentReducedCost[q]);
        Enorm += exp(kappaCol*currentPotentialEnergy[q]);
    }
    assert(currentPotentialEnergy == potentialEnergies.front());
    assert(Enorm == energyNormalisation.front());
    double cumulativeP = 0.0;
    std::vector<double> cdist(simplex.nNonBasic() + 1);
    cdist[0] = 0.0;
    for(int q=1; q <= simplex.nNonBasic(); ++q) {
        cumulativeP += exp(kappaCol* currentPotentialEnergy[q]);
        cdist[q] = cumulativeP;
    }
//    std::cout << "cdf = " << cdf << std::endl;
//    std::cout << "cdist = " << cdist << std::endl;
    assert(cdist == cdf);
}

// returns -1, 0 or +1
//int NormalisedPotentialEnergyPivot::sign(double x) {
//    return (x > tol) - (x < -tol);
//}

//double NormalisedPotentialEnergyPivot::potentialEnergy(int j, const std::vector<double> &reducedCost) {
//    return potentialEnergy(simplex.isAtUpperBound(j), reducedCost[j]);
//}
