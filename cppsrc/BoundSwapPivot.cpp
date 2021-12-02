//
// Created by daniel on 27/11/2021.
//

#include "BoundSwapPivot.h"
#include "SimplexMCMC.h"

BoundSwapPivot::BoundSwapPivot(SimplexMCMC &simplex)
        :
        ProposalPivot(simplex,0,0),
//kappaCol(std::max(std::log(p1*simplex.nNonBasic()),1.0)),  // kappaCol(3.5) {
//kappaCol(std::max(std::log((simplex.nNonBasic()-1.0)/(1.0/p1 - 1.0)),1.0)),  // kappaCol(3.5) {
        cdf(simplex.nNonBasic() + 1)
{ }


// called from SimplexMCMC
void BoundSwapPivot::init() {
    initBasis();
    randomiseBounds();
    calculateTableau();
    col.resize(simplex.nBasic()+1);
    deltaj = 0.0;
    currentInfeasibilityCosts.resize(simplex.nBasic()+1, 0.0);
    currentInfeasibilityCosts[0] = 0.0;
    currentPotentialEnergies.resize(simplex.nNonBasic()+1, 0.0);
    currentPotentialEnergies[0] = 0.0;
    recalculateInfeasibilityCost();
    recalculateReducedCosts();
    recalculatePotentials();
    cdf[0] = 0.0;
    recalculateCDF();
}


// TODO: generalise this
void BoundSwapPivot::initBasis() {
//    std::cout << "Number of fixed vars = " << simplex.kProbTokSim.size() - simplex.kSimTokProb.size() << std::endl;
//    int nAux = 0;
//    for(int i = 1; i< simplex.nBasic(); ++i) {
//        if(simplex.isAuxiliary(simplex.head[i])) {
//            ++nAux;
//        } else {
//            std::cout << "Pivoted-in Event = " << Event<PredPreyAgent<4>>(simplex.kSimTokProb[simplex.head[i]]-simplex.nBasic()) << std::endl;
//        }
//    }
//    std::cout << "Number of aux = " << nAux << std::endl;
//    std::cout << "Number of basic = " << simplex.nBasic() << std::endl;
    for(int j = simplex.nNonBasic(); j>0; --j) {
        if(isInPredPreyPreferredBasis(simplex.head[j+simplex.nBasic()])) {
            std::vector<double> col = simplex.tableauCol(j);
            int i=0;
            while(i <= simplex.nBasic() && (col[i]==0.0 || simplex.isAuxiliary(simplex.head[i]))) ++i;
            if(i <= simplex.nBasic()) {
//                std::cout << "Pivoting in " << i << " " << j << std::endl;
                simplex.pivot(i, j, col, simplex.beta[i] == 1.0);
            }
        }
    }
}


// For debug only: set bounds of nonBasics to random vals
void BoundSwapPivot::randomiseBounds() {
    for(int j=1; j<simplex.nNonBasic(); ++j) {
        simplex.flag[j] = Random::nextBool();
    }
    simplex.evalBeta();
}

void BoundSwapPivot::calculateTableau() {
    tableauCols.reserve(simplex.nNonBasic()+1);
    tableauCols.push_back(glp::SparseVec());
    for(j=1; j <= simplex.nNonBasic(); ++j) {
        tableauCols.push_back(glp::SparseVec(simplex.tableauCol(j)));
    }
}

const BoundSwapPivot &BoundSwapPivot::nextProposal() {
    // update cache
    if(simplex.lastSampleWasAccepted) {
        // remove last state's cache
        bool infeasibilityCostHasChanged = false;
        if(deltaj != 0.0) {
            infeasibilityCostHasChanged = recalculateInfeasibilityCost();
        }

        // need to recalculate potential energy if
        //   - infeasibility cost has changed
        //   - we did a real pivot on a column with non-zero reduced cost
        if (infeasibilityCostHasChanged || (i > 0 && currentReducedCosts[j] != 0.0)) {
            recalculateReducedCosts();
            recalculatePotentials();
            recalculateCDF();
        }
    }
//    checkCurrentCacheIsValid();
    chooseCol();
    chooseRow();
    calcAcceptanceContrib();
    return *this;
}

bool BoundSwapPivot::recalculateInfeasibilityCost() {
    bool hasChanged = false;
    for(int i=1; i<=simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        double infeasibility = infeasibilityGradient(simplex.beta[i], simplex.l[k], simplex.u[k]);
        if(infeasibility != currentInfeasibilityCosts[i]) {
            currentInfeasibilityCosts[i] = infeasibility;
            hasChanged = true;
        }
    }
    return hasChanged;
}


void BoundSwapPivot::recalculateReducedCosts() {
    std::vector<double> CBI(currentInfeasibilityCosts); // will eventually be proposed C'B'^-1
    simplex.btran(CBI); // in-place post-multiplication by B^-1
    currentReducedCosts = simplex.piTimesMinusN(CBI);
};

void BoundSwapPivot::recalculatePotentials() {
    for (int q = 1; q <= simplex.nNonBasic(); ++q) {
        currentPotentialEnergies[q] = potentialEnergy(simplex.isAtUpperBound(q), currentReducedCosts[q]);
    }
};

void BoundSwapPivot::recalculateCDF() {
    double cumulativeP = 0.0;
    for(int q=1; q <= simplex.nNonBasic(); ++q) {
        cumulativeP += exp(kappaCol* currentPotentialEnergies[q]);
        cdf[q] = cumulativeP;
    }
}



bool BoundSwapPivot::isInPredPreyPreferredBasis(int k) {
    int kProb = simplex.kSimTokProb[k];
    return kProb>simplex.nBasic() && (kProb-simplex.nBasic()-1)%PredPreyAgentBase::actDomainSize() == PredPreyAgentBase::DIE;
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
void BoundSwapPivot::chooseCol() {
    auto chosenJ = std::lower_bound(
            cdf.begin(),
            cdf.end(),
            Random::nextDouble(0.0, cdf[simplex.nNonBasic()])
    );
//    setCol(chosenJ - cdf.begin());
    j = chosenJ - cdf.begin();
    clearCol();
    nonZeroRows = tableauCols[j].indices;
    for(int nzz =0; nzz<nonZeroRows.size(); ++nzz) {
        int nzi = nonZeroRows[nzz];
        col[nzi] = tableauCols[j].values[nzz];
    }
}



// Swap with probaability
// exp(kappaCol*destinationPotential + kappaRow*destinationInfeasibility)/
// (exp(kappaCol*currentPotential + kappaRow*cuurentInfeasibility) + exp(kappaCol*destinationPotential + kappaRow*destinationInfeasibility))
// =
// exp(kappaCol*DeltaPotential + kappaRow*DeltaInfeasibility)/(1 + exp(kappaCol*DeltaPotential + kappaRow*DeltaInfeasibility))
void BoundSwapPivot::chooseRow() {
    double dj = simplex.isAtUpperBound(j)?-1.0:1.0; // Assumes Fermionic
    double deltaInfeasibility = 0.0;
    double swapReducedCost = 0.0;
    for(int nzz = 0; nzz < tableauCols[j].sparseSize(); ++nzz) {
        int nzi = tableauCols[j].indices[nzz];
        int rowk = simplex.head[nzi];
        double bi = simplex.beta[nzi];
        double upperBound = simplex.u[rowk];
        double lowerBound = simplex.l[rowk];
        double Mij = tableauCols[j].values[nzz];
        double swapbi = bi + Mij * dj;
        deltaInfeasibility += infeasibility(swapbi,lowerBound, upperBound) - infeasibility(bi, lowerBound, upperBound);
        swapReducedCost += Mij*infeasibilityCostFn(swapbi,lowerBound, upperBound);
    }
//    assert(swapReducedCost == colInfeasibilityGradient(dj));
    double swapPotential = potentialEnergy(!simplex.isAtUpperBound(j), swapReducedCost);
    double deltaPotential = swapPotential - currentPotentialEnergies[j];
    double deltaProb = exp(kappaCol*deltaPotential + kappaRow*deltaInfeasibility);

    bool doSwap = (Random::nextDouble()*(1.0 + deltaProb) < deltaProb);

    leavingVarToUpperBound = simplex.isAtUpperBound(j) xor doSwap;
    i = -1; // column does bound swap.
    deltaj = doSwap?dj:0.0;
}


void BoundSwapPivot::applyProposal() {

}



// Calculates log acceptance contribution which is kappaCol times
// the proposed change in the potential energy, ignoring column j.
void BoundSwapPivot::calcAcceptanceContrib() {
    logAcceptanceContribution = 0.0;
}


