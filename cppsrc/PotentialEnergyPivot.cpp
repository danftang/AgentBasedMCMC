//
// Created by daniel on 06/07/2021.
//

#include <cfloat>

#include "Random.h"
#include "StlStream.h"
#include "PotentialEnergyPivot.h"
#include "SimplexMCMC.h"

PotentialEnergyPivot::PotentialEnergyPivot(SimplexMCMC &simplex)
:
ProposalPivot(simplex,0,0),
kappaCol(std::max(std::log(p1*simplex.nNonBasic()),1.0)),  // kappaCol(3.5) {
cdf(simplex.nNonBasic() + 1)
{
    // setup initial cached values
    initCache();
}

void PotentialEnergyPivot::initCache() {
    infeasibilityCosts.push_back(infeasibilityCost());
    calcProposedReducedCosts();
    calcProposedEnergies();
    cdf[0] = 0.0;
    recalculateCDF();
}

const PotentialEnergyPivot &PotentialEnergyPivot::nextProposal() {
    // reset cache
    if(simplex.lastSampleWasAccepted) {
        if(infeasibilityCosts.size() > 1)   infeasibilityCosts.pop_front();
        if(reducedCosts.size() > 1) reducedCosts.pop_front();
        if(potentialEnergies.size() > 1) {
            potentialEnergies.pop_front();
            recalculateCDF();
        }
        if(totalPotentials.size() >1) totalPotentials.pop_front();
    } else {
        if(infeasibilityCosts.size() > 1)   infeasibilityCosts.pop_back();
        if(reducedCosts.size() > 1) reducedCosts.pop_back();
        if(potentialEnergies.size() > 1)    potentialEnergies.pop_back();
        if(totalPotentials.size() >1) totalPotentials.pop_back();
    }
//    checkCurrentCacheIsValid();
    chooseCol();
    chooseRow();
    calcAcceptanceContrib();
    return *this;
}


void PotentialEnergyPivot::recalculateCDF() {
    double cumulativeP = 0.0;
    const std::vector<double> &currentPotentials = potentialEnergies.front();
    for(int q=1; q <= simplex.nNonBasic(); ++q) {
//        double colPotential = potentialEnergy(q, currentReducedCost);
//        Ep += colPotential;
        cumulativeP += exp(kappaCol* currentPotentials[q]);
        cdf[q] = cumulativeP;
    }
//    assert(cumulativeP > simplex.nNonBasic()); // there should be at least one high energy column
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
void PotentialEnergyPivot::chooseCol() {
    auto chosenJ = std::lower_bound(
            cdf.begin(),
            cdf.end(),
            Random::nextDouble(0.0, cdf[simplex.nNonBasic()])
        );
    setCol(chosenJ - cdf.begin());
}



void PotentialEnergyPivot::chooseRow() {
    std::multimap<double, int> transitions = getPivotsByDeltaJ(); // from delta_j to LogPMF-index.


    // now populate pivotPMF by going from lowest to highest delta_j in order
    std::vector<double> pivotPMF(nonZeroRows.size() * 2 + 2, DBL_MAX); // index is (2*nonZeroRowIndex + toUpperBound), final two are lower and upper bound swap, value is probability mass
    double lastDj = transitions.begin()->first;
    double infeas = 0.0; //infeasibility(lastDj) - infeasibility(0.0);
    double infeasMin = DBL_MAX;
    double dDf_dDj = colInfeasibilityGradient(lastDj - 2.0*tol);
    for(auto [Dj, pmfIndex] : transitions) {
        infeas += (Dj - lastDj) * dDf_dDj;
        lastDj = Dj;
        if(isActive(pmfIndex)) {
            pivotPMF[pmfIndex] = infeas;
            if(infeas < infeasMin) infeasMin = infeas;
        }
//        std::cout << "deltaj = " << Dj << " pivotProb = " << pivotPMF[pmfIndex] << " infeas = " << infeas << " dDf_dDj = " << dDf_dDj << " Feasibility = " << infeasibility(Dj) << std::endl;
        if(pmfIndex < nonZeroRows.size()*2) {
            dDf_dDj += fabs(col[nonZeroRows[pmfIndex / 2]]);
        } else {
            dDf_dDj += 1.0;
        }
    }
    assert(infeasMin != DBL_MAX);

    // rescale and take exponential
    for(int i=0; i<pivotPMF.size(); ++i) {
        pivotPMF[i] = exp(kappaRow*(pivotPMF[i] - infeasMin));
    }
//    pivotPMF[2*nonZeroRows.size() + simplex.isAtUpperBound(j)] *= 0.01; // preference against null pivot if other possibilities exist

    // nextIntFromDiscrete row

    setToPivotIndex(Random::nextIntFromDiscrete(pivotPMF.begin(), pivotPMF.end()));

    //    std::cout << "Chose pivot with prob " << pivotPMF[pivotChoice] << " Delta = " << deltaj << std::endl;
//  std::cout << "Proposing pivot from " << !sourceObjectiveIsZero << " to " << destinationObjective << std::endl;
}



// Calculates log acceptance contribution which is kappaCol times
// the proposed change in the potential energy, ignoring column j.
//
// To calculate the post-pivot potential, we need to multiply the post-pivot
// cost C' by the post-pivot basis inverse, B'^-1, then multiply by N.
// We can express C'B'^-1 as (C'E)B^-1 where B^-1 is the pre-pivot basis inverse
// and E is the elementary transform matrix EB^-1 = B'^-1.
//
//
void PotentialEnergyPivot::calcAcceptanceContrib() {

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

    // logAcceptanceContribution is now -k times the change in total potential energy in all cols q!=j
    if(totalPotentials.size() > 1) {
        logAcceptanceContribution = kappaCol * (
                (totalPotentials.front() - potentialEnergies.front()[j]) -
                (totalPotentials.back()  - potentialEnergies.back()[j])
                );
    } else {
        logAcceptanceContribution = 0.0;
    }
//    std::cout << "Log acceptance contribution = " << logAcceptanceContribution << std::endl;
}


void PotentialEnergyPivot::calcProposedInfeasibilityCosts() {
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
        proposedInfeasibilityCost[p] = infeasibilityGradient(postPivotBetap, simplex.l[k], simplex.u[k]);
        if(currentInfeasibilityCost[p] != proposedInfeasibilityCost[p]) hasChanged = true;
    }
    if(hasChanged) infeasibilityCosts.push_back(std::move(proposedInfeasibilityCost));
}


void PotentialEnergyPivot::calcProposedReducedCosts() {
    std::vector<double> proposedCBI(infeasibilityCosts.back()); // will eventually be proposed C'B'^-1
    // if we're actually pivoting (i>0) update infeasibilityCost:
    // proposedCBI = proposedCBI * E, where E is the elementary transform matrix (see above).
    if (i > 0) {
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
        reducedCosts.push_back(std::move(proposedReducedCost));
    } else {
        simplex.btran(proposedCBI); // in-place post-multiplication by B^-1
        std::vector<double> proposedReducedCost = simplex.piTimesMinusN(proposedCBI);
        reducedCosts.push_back(std::move(proposedReducedCost));
    }
}

void PotentialEnergyPivot::calcProposedEnergies() {
    double proposedEp = 0.0;
    std::vector<double> proposedPotential(simplex.nNonBasic()+1);
    std::vector<double> &proposedReducedCost = reducedCosts.back();
    for (int q = 1; q <= simplex.nNonBasic(); ++q) {
        bool proposedQIsAtUpperBound = (q == j)?leavingVarToUpperBound:simplex.isAtUpperBound(q);
        proposedEp += proposedPotential[q] = potentialEnergy(proposedQIsAtUpperBound, proposedReducedCost[q]);
    }
    potentialEnergies.push_back(std::move(proposedPotential));
    totalPotentials.push_back(proposedEp);
}

// Sanity check for calcAcceptanceContrib()
void PotentialEnergyPivot::calcAcceptanceContribCheck() {

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
    totalPotentials.push_back(postPivotEp);
    std::vector<double> revCol = simplex.tableauCol(j);
    simplex.pivot(i, j, revCol, jIsAtUpperBound);

    logAcceptanceContribution = kappaCol*(
            (totalPotentials.front() - potentialEnergies.front()[j]) -
            (totalPotentials.back() - potentialEnergies.back()[j]));
}

void PotentialEnergyPivot::checkCurrentCacheIsValid() {
    std::vector<double> currentCost = infeasibilityCost();
    assert(currentCost == infeasibilityCosts.front());
    simplex.btran(currentCost);
    std::vector<double> currentReducedCost = simplex.piTimesMinusN(currentCost);
    assert(currentReducedCost == reducedCosts.front());
    double Ep = 0.0;
    std::vector<double> currentPotentialEnergy(simplex.nNonBasic()+1);
    currentPotentialEnergy[0] = 0.0;
    for(int q=1; q <= simplex.nNonBasic(); ++q) {
        currentPotentialEnergy[q] = potentialEnergy(simplex.isAtUpperBound(q),currentReducedCost[q]);
        Ep += currentPotentialEnergy[q];
    }
    assert(currentPotentialEnergy == potentialEnergies.front());
    assert(Ep == totalPotentials.front());
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
//int PotentialEnergyPivot::sign(double x) {
//    return (x > tol) - (x < -tol);
//}

//double PotentialEnergyPivot::potentialEnergy(int j, const std::vector<double> &reducedCost) {
//    return potentialEnergy(simplex.isAtUpperBound(j), reducedCost[j]);
//}
