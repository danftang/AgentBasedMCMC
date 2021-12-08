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
        currentE(simplex.nBasic()+1),
        currentDeltaE(simplex.nNonBasic()+1),
        feasibleEnergy(simplex.nVars()+1),
        colPMF(simplex.nNonBasic() + 1)
{ }


// called from SimplexMCMC
void BoundSwapPivot::init() {
    initBasis();
//    randomiseBounds();
    calculateTableau();
    deltaj = 0.0;

    // calculate feasible energies
    for(int kSim=1; kSim<=simplex.nVars(); ++kSim) {
        int kProb = simplex.kSimTokProb[kSim];
        if(kProb > simplex.nBasic()) {
            auto event = Event<PredPreyAgent<8>>(kProb-simplex.nBasic()); // TODO: make this generic
            feasibleEnergy[kSim] = log(event.agent().marginalTimestep()[event.act()])/kappa;
        } else {
            feasibleEnergy[kSim] = 0.0;
        }
    }

    // calculate row energies
    for(int i=1; i<=simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        currentE[i] = energy(i, simplex.beta[i]);
    }

    currentDeltaE = calcDeltaE();

    // calculate column PMF
    std::vector<double> colDistribution(simplex.nNonBasic()+1);
    colDistribution[0] = 0.0;
    for(int j=1; j<=simplex.nNonBasic(); ++j) {
        colDistribution[j] = currentDeltaE[j]>0.0?exp(kappa*currentDeltaE[j]):1.0;
    }
    colPMF.setAll(colDistribution);


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
            while(i <= simplex.nBasic() && (col[i]==0.0 || simplex.isAuxiliary(simplex.head[i]) ||
                    isInPredPreyPreferredBasis(simplex.head[i]))) ++i;
            if(i <= simplex.nBasic()) {
//                std::cout << "Pivoting in preferred basis " << i << " " << j << std::endl;
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
    chooseBound();
    if(deltaj != 0.0) {
//        std::cout << j << " " << deltaj << " " << deltaFeasibleEnergy() << " "
//                  << (currentDeltaE[j] - deltaInfeasibility()) << std::endl;
//        assert(fabs(deltaFeasibleEnergy() - (currentDeltaE[j] - deltaInfeasibility())) < 1e-8);
        logAcceptanceContribution = -kappa * (currentDeltaE[j] - deltaInfeasibility());
    } else {
        logAcceptanceContribution = 0.0;
    }
//    std::cout << "Proposing " << j << " " << deltaj << " " << leavingVarToUpperBound << std::endl;
    return *this;
}

//bool BoundSwapPivot::recalculateInfeasibilityCost() {
//    bool hasChanged = false;
//    for(int i=1; i<=simplex.nBasic(); ++i) {
//        double infeasibility = infeasibilityCostFn(i);
//        if(infeasibility != currentInfeasibilityCosts[i]) {
//            currentInfeasibilityCosts[i] = infeasibility;
//            hasChanged = true;
//        }
//    }
//    return hasChanged;
//}

// assuming the last proposal was accepted, beta has been
// updated accordingly and changedDeltaE constains the changes
// to deltaE caused by the pivot.
void BoundSwapPivot::updateAllCosts() {
    const glp::SparseVec &col = tableauCols[j];
    for(int nzi=0; nzi<col.sparseSize(); ++nzi) {
        int i = col.indices[nzi];
        currentE[i] = energy(i, simplex.beta[i]);
    }
    for(const auto [changedj, newDeltaj]: changedDeltaE) {
        currentDeltaE[changedj] = newDeltaj;
        colPMF[changedj] = newDeltaj>0.0?exp(kappa*newDeltaj):1.0;
    }
}


// calculates the changes to deltaE by col if we were to pivot on column j;
void BoundSwapPivot::calcDeltaEChanges() {
    changedDeltaE.clear();
    const glp::SparseVec &col = tableauCols[j];
    changedDeltaE[j] = -currentDeltaE[j];
    double dj = simplex.isAtUpperBound(j)?-1.0:1.0;
    for(int nzi=0; nzi<col.sparseSize(); ++nzi) {
        int i = col.indices[nzi]; // row that is affected by the change
        double oldbetai = simplex.beta[i];
        double newbetai = oldbetai + col.values[nzi]*dj; // the value of beta_i after swapping col j;
        double newEi = energy(i, newbetai); // new Energy of row i
        const glp::SparseVec &row = tableauRows[i];
        double deltaEi = newEi - currentE[i];
        for (int nzi = 0; nzi < row.sparseSize(); ++nzi) {
            int updatej = row.indices[nzi];
            if(updatej != j) {
                changedDeltaE.try_emplace(updatej, currentDeltaE[updatej]);
                double deltabetaiswapj = row.values[nzi] * (simplex.isAtUpperBound(updatej) ? -1.0
                                                                                            : 1.0); // the change in beta_i on swapping col updatej
                changedDeltaE[updatej] +=
                        (energy(i, newbetai + deltabetaiswapj) - newEi)
                        -(energy(i, oldbetai + deltabetaiswapj) - currentE[i]);
            }
        }
    }
}

// Calculates changes in energy by row, assumes currentE is valid.
std::vector<double> BoundSwapPivot::calcDeltaE() {
    std::vector<double> deltaE(simplex.nNonBasic()+1);
    deltaE[0] = 0.0;
    for(int j=1; j<=simplex.nNonBasic(); ++j) {
        double dj = simplex.isAtUpperBound(j)?-1.0:1.0;
        deltaE[j] = dj*feasibleEnergy[simplex.head[j+simplex.nBasic()]];
        const glp::SparseVec &col = tableauCols[j];
        for(int nzi = 0; nzi<col.sparseSize(); ++nzi) {
            int i = col.indices[nzi];
            double Mij = col.values[nzi];
            double oldbetai = simplex.beta[i];
            double newbetai = oldbetai + col.values[nzi]*dj; // the value of beta_i after swapping col j;
            double newEi = energy(i, newbetai); // new Energy of row i
            deltaE[j] += newEi - currentE[i];
        }
    }
    return deltaE;
}


//void BoundSwapPivot::recalculateCDF() {
//    double cumulativeP = 0.0;
//    std::vector<double> P(simplex.nNonBasic()+1);
//    P[0] = 0.0;
//    for(int q=1; q <= simplex.nNonBasic(); ++q) {
//        P[q] = exp(kappaCol* potentialEnergy(simplex.isAtUpperBound(q), currentReducedCosts[q]));
//    }
//    colPMF.setAll(P);
//}



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
    j = colPMF(Random::gen);
}


// Swap with probaability
// exp(kappaRow*destinationInfeasibility)P(chooseCol|destination)/
// (exp(kappaRow*cuurentInfeasibility)P(chooseCol|currentState) + exp(kappaRow*destinationInfeasibility)P(chooseCol|destination))
// =
// exp(kappaCol*DeltaPotential + kappaRow*DeltaInfeasibility)/(1 + exp(kappaCol*DeltaPotential + kappaRow*DeltaInfeasibility))
void BoundSwapPivot::chooseBound() {
    calcDeltaEChanges();

    double changeInSum = 0.0;
    for(const auto [changedj, newDeltaj]: changedDeltaE) {
        double oldPj = colPMF[changedj];
        changeInSum += (newDeltaj>0.0?exp(kappa*newDeltaj):1.0) - oldPj;
    }
    double ratioOfSums = colPMF.sum()/(colPMF.sum() + changeInSum);
    bool doSwap;
    if(ratioOfSums >= 1.0) {
        doSwap = true;
    } else {
        doSwap = (Random::nextDouble() < ratioOfSums);
    }
    if(doSwap) {
        leavingVarToUpperBound = !simplex.isAtUpperBound(j);
        deltaj = simplex.isAtUpperBound(j)?-1.0:1.0;
    } else {
        leavingVarToUpperBound = simplex.isAtUpperBound(j);
        deltaj = 0.0;
    }
}


double BoundSwapPivot::energy(int i, double b) {
    int k = simplex.head[i];
    return energy(b, feasibleEnergy[k], simplex.l[k], simplex.u[k]);
}



//double BoundSwapPivot::infeasibilityCostFn(int i) {
//    int k = simplex.head[i];
//    return infeasibilityCostFn(simplex.beta[i], simplex.l[k], simplex.u[k]);
//}


// Calculates the change in feasible energy given the current
// values j and deltaj
double BoundSwapPivot::deltaFeasibleEnergy() {
    double dEf = deltaj*feasibleEnergy[simplex.head[j+simplex.nBasic()]];
    const glp::SparseVec &col = tableauCols[j];
    for(int nzi=0; nzi<col.sparseSize(); ++nzi) {
        int i = col.indices[nzi];
        double oldbetai = simplex.beta[i];
        double newbetai = oldbetai + col.values[nzi]*deltaj;
        if(oldbetai > 1.0) oldbetai = 1.0;
        if(oldbetai < 0.0) oldbetai = 0.0;
        if(newbetai > 1.0) newbetai = 1.0;
        if(newbetai < 0.0) newbetai = 0.0;
        dEf += (newbetai - oldbetai)*feasibleEnergy[simplex.head[i]];
    }
    return dEf;
}


double BoundSwapPivot::deltaInfeasibility() {
    double dI = 0.0;
    const glp::SparseVec &col = tableauCols[j];
    for(int nzi=0; nzi<col.sparseSize(); ++nzi) {
        int i = col.indices[nzi];
        int k = simplex.head[i];
        double upperBound = simplex.u[k];
        double lowerBound = simplex.l[k];
        double oldbetai = simplex.beta[i];
        double newbetai = oldbetai + col.values[nzi]*deltaj;
        dI += infeasibility(newbetai, lowerBound, upperBound) - infeasibility(oldbetai, lowerBound, upperBound);
    }
    return dI;
}


void BoundSwapPivot::checkCosts() {
    for(int i=1; i<=simplex.nBasic();++i) {
        assert(fabs(currentE[i] - energy(i, simplex.beta[i])) < 1e-8);
    }

    std::vector<double> deltaE = calcDeltaE();
    for(int j=1; j<=simplex.nNonBasic(); ++j) {
        if(fabs(deltaE[j] - currentDeltaE[j]) > 1e-8) {
            std::cout << j << " " << simplex.nNonBasic() << " " << deltaE[j] << " " << currentE[j] << std::endl;
            assert(false);
        }
        double Pj = deltaE[j]>0.0?exp(kappa*deltaE[j]):1.0;
        if(fabs(colPMF[j] - Pj) > 1e-8) {
            std::cout << j << " " << simplex.nNonBasic() << " " << (double)colPMF[j] << " " << Pj << std::endl;
            assert(false);
        }
    }
}
