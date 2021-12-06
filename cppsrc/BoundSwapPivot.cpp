//
// Created by daniel on 27/11/2021.
//

#include "BoundSwapPivot.h"
#include "SimplexMCMC.h"

BoundSwapPivot::BoundSwapPivot(SimplexMCMC &simplex)
        :
        simplex(simplex),
        i(-1),
        j(0),
        logAcceptanceContribution(0.0),
        tableauRows(simplex.nBasic()+1,glp::SparseVec()),
        cdf(simplex.nNonBasic() + 1)
{ }


// called from SimplexMCMC
void BoundSwapPivot::init() {
    initBasis();
//    randomiseBounds();
    calculateTableau();
    deltaj = 0.0;
    currentInfeasibilityCosts.resize(simplex.nBasic()+1, 0.0);
    currentInfeasibilityCosts[0] = 0.0;
    recalculateInfeasibilityCost();
    currentReducedCosts = reducedCosts();
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
    for(int j=1; j <= simplex.nNonBasic(); ++j) {
        tableauCols.push_back(glp::SparseVec(simplex.tableauCol(j)));
        const glp::SparseVec &col = tableauCols[j];
        for(int nzi=0; nzi<col.sparseSize(); ++nzi) tableauRows[col.indices[nzi]].insert(j, col.values[nzi]);
    }
}

BoundSwapPivot &BoundSwapPivot::nextProposal() {
    if(simplex.lastSampleWasAccepted && deltaj != 0.0) updateAllCosts();
//    checkCosts();
    chooseCol();
    chooseRow();
    return *this;
}

bool BoundSwapPivot::recalculateInfeasibilityCost() {
    bool hasChanged = false;
    for(int i=1; i<=simplex.nBasic(); ++i) {
        double infeasibility = simplex.infeasibilityGradient(i);
        if(infeasibility != currentInfeasibilityCosts[i]) {
            currentInfeasibilityCosts[i] = infeasibility;
            hasChanged = true;
        }
    }
    return hasChanged;
}


void BoundSwapPivot::updateAllCosts() {
    const glp::SparseVec &col = tableauCol();
    std::set<int> changedCols;
    changedCols.insert(j);
    for(int nzi=0; nzi<col.sparseSize(); ++nzi) {
        int i = col.indices[nzi];
        double infeasibility = infeasibilityCostFn(i);
        if(infeasibility != currentInfeasibilityCosts[i]) {
            // update reduced costs
            double deltaInfi = infeasibility - currentInfeasibilityCosts[i];
            const glp::SparseVec &row = tableauRows[i];
            for(int nzi=0; nzi < row.sparseSize(); ++nzi) {
                int j = row.indices[nzi];
                currentReducedCosts[j] += deltaInfi * row.values[nzi]; // not worried about loss of precision as always dealing with integers
                changedCols.insert(j);
            }
            currentInfeasibilityCosts[i] = infeasibility;
        }
    }
    for(int j : changedCols) {
        cdf[j] = exp(kappaCol*potentialEnergy(simplex.isAtUpperBound(j), currentReducedCosts[j]));
    }
}


std::vector<double> BoundSwapPivot::reducedCosts() {
    std::vector<double> CBI(currentInfeasibilityCosts); // will eventually be proposed C'B'^-1
    simplex.btran(CBI); // in-place post-multiplication by B^-1
    return simplex.piTimesMinusN(CBI);
};


void BoundSwapPivot::recalculateCDF() {
    double cumulativeP = 0.0;
    std::vector<double> P(simplex.nNonBasic()+1);
    P[0] = 0.0;
    for(int q=1; q <= simplex.nNonBasic(); ++q) {
        P[q] = exp(kappaCol* potentialEnergy(simplex.isAtUpperBound(q), currentReducedCosts[q]));
    }
    cdf.setAll(P);
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
    j = cdf(Random::gen);
}


// Swap with probaability
// exp(kappaRow*destinationInfeasibility)P(chooseCol|destination)/
// (exp(kappaRow*cuurentInfeasibility)P(chooseCol|currentState) + exp(kappaRow*destinationInfeasibility)P(chooseCol|destination))
// =
// exp(kappaCol*DeltaPotential + kappaRow*DeltaInfeasibility)/(1 + exp(kappaCol*DeltaPotential + kappaRow*DeltaInfeasibility))
void BoundSwapPivot::chooseRow() {
    double dj = simplex.isAtUpperBound(j)?-1.0:1.0; // Assumes Fermionic
    double deltaInfeasibility = 0.0;
    double swapReducedCost = 0.0;
    for(int nzz = 0; nzz < tableauCols[j].sparseSize(); ++nzz) {
        int nzi = tableauCols[j].indices[nzz];
        double Mij = tableauCols[j].values[nzz];
        int rowk = simplex.head[nzi];
        double bi = simplex.beta[nzi];
        double upperBound = simplex.u[rowk];
        double lowerBound = simplex.l[rowk];
        double swapbi = bi + Mij * dj;
        deltaInfeasibility += infeasibility(swapbi,lowerBound, upperBound) - infeasibility(bi, lowerBound, upperBound);
        swapReducedCost += Mij*infeasibilityCostFn(swapbi,lowerBound, upperBound);
    }
//    assert(swapReducedCost == colInfeasibilityGradient(dj));
    double swapPotential = potentialEnergy(!simplex.isAtUpperBound(j), swapReducedCost);
    double swapChooseColScore = exp(kappaCol*swapPotential);
    double currentChooseColScore = cdf[j];
//    double pChooseColRatio = swapChooseColScore*cdf.sum()/(currentChooseColScore*(cdf.sum() + swapChooseColScore - currentChooseColScore));
    double pChooseColRatio = swapChooseColScore/currentChooseColScore;
    double deltaProb = exp(kappaRow*deltaInfeasibility)*pChooseColRatio;

    bool doSwap = (Random::nextDouble()*(1.0 + deltaProb) < deltaProb);

//    std::cout << doSwap << " ChooseColRatio = " << pChooseColRatio << " " << deltaProb/(1.0 + deltaProb) << std::endl;

    leavingVarToUpperBound = simplex.isAtUpperBound(j) xor doSwap;
//    i = -1; // column does bound swap.
    deltaj = doSwap?dj:0.0;
}



double BoundSwapPivot::infeasibilityCostFn(int i) {
    int k = simplex.head[i];
    return infeasibilityCostFn(simplex.beta[i], simplex.l[k], simplex.u[k]);
}

void BoundSwapPivot::checkCosts() {
    for(int i=1; i<=simplex.nBasic();++i) {
        assert(fabs(currentInfeasibilityCosts[i] - infeasibilityCostFn(i)) < 1e-8);
    }

    std::vector<double> rc = reducedCosts();
    for(int j=1; j<=simplex.nNonBasic(); ++j) {
        assert(fabs(rc[j] - currentReducedCosts[j]) < 1e-8);
        if(fabs(cdf[j] - exp(kappaCol*potentialEnergy(simplex.isAtUpperBound(j), rc[j]))) > 1e-8) {
            std::cout << j << " " << simplex.nNonBasic() << " " << (double)cdf[j] << " " << exp(kappaCol*potentialEnergy(simplex.isAtUpperBound(j), rc[j])) << " " << fabs(cdf[j] - exp(kappaCol*potentialEnergy(simplex.isAtUpperBound(j), rc[j]))) << std::endl;
            assert(false);
        }
    }
}
