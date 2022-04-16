// Samples from a lattice of points in a plane represented as a set of basis vectors, stored in sparse format.
// The valid points of the lattice are contained within a hyperrectangle with one corner at the origin
// and a diagonally opposite corner at "H", so that the valid points can be defined as the set
//
// { A : A = NX_N + F, L <= A <= U }
//
// where the columns of N are the basis vectors, and the elements of $X_N$ are integers.
//
// The probability of each point must be expressible as a product of unnormalised marginals
//
// P(A) = K\prod_i P_i(A_i) \prod_i I(A_i)
//
// where I(X_i) = { 1               if L_i <= X_i <= U_i,
//                { e^{k(L_i-X_i}   if X_i < L_i
//                { e^{k(H_i-Xi)}   if X_i > H_i
//
// Created by daniel on 07/03/2022.
//

#ifndef ABMCMC_SPARSEBASISSAMPLER_H
#define ABMCMC_SPARSEBASISSAMPLER_H

#include <utility>
#include <vector>
#include "TableauNormMinimiser.h"
#include "MutableCategoricalArray.h"
#include "MCMCStatistics.h"
#include "Random.h"

// TODO: need to generalise marginalGradE to general function, get rid of X (now in importance function),
//  integrate importance function and finish construction from WeightedFactoredConvexDistribution.

template<class T>
class SparseBasisSampler {
public:
//    static constexpr double   kappa = 2.0; // (+ve) rate of decay of probability outside the hyper-rectangle

    double   kappa; // (+ve) rate of decay of probability outside the hyper-rectangle


    // tableau
    std::vector<SparseVec<T>>   cols;   // the basis vectors
    std::vector<SparseVec<T>>   rows;   //
    std::vector<T>              L;      // lower bounds by row
    std::vector<T>              U;      // upper bounds by row
    std::vector<std::function<double(T)>> factors;

    std::vector<T>              X;      // the current point on the lattice by row.
    std::vector<T>              delta;  // direction of perturbation of each basis by column (assume act Fermionicity)
    MutableCategoricalArray     basisDistribution; // distribution of proposing to update basis by column
    std::vector<double>         currentDE;      // change in (-ve) energy by flipping the value of the j'th column
    std::vector<double>         currentE;       // current (-ve) energy of i'th row
    T                           currentInfeasibility;
    double                      currentLogImportance; // real prob over estimated prob of current state (=1 if infeasible)

    std::unique_ptr<PerturbableFunction<T,double>>    importanceFunc;

    MCMCStatistics              stats;

    class Proposal {
    public:
        Proposal(SparseBasisSampler<T> &parent);

        int proposedj;                           // column index of proposed basis
        const std::vector<int> &changedXIndices; // indices of X that are changed in the proposal
        std::vector<T>          changedXValues;  // new values of elements of X, by sparse index of j'th basis
        std::map<int,double>    changedDE;       // change in Delta-Energy by column
        T                       infeasibility;   // proposed infeasibility

        double calcRatioOfSums(const MutableCategoricalArray &basisDistribution);
    };

    SparseBasisSampler(const TableauNormMinimiser<T> &tableau, const std::vector<std::function<double(T)>> &tableauFactors, std::unique_ptr<PerturbableFunction<T,double>> &&importanceFunc, const std::vector<T> &initialState, double Kappa);

    SparseBasisSampler(const WeightedFactoredConvexDistribution<T> &distribution, const std::vector<T> & initialState, double Kappa);

    SparseBasisSampler(const WeightedFactoredConvexDistribution<T> &distribution, double Kappa);

    const std::vector<T> &operator()();     // take the next sample

    void findInitialFeasibleSolution();     // make the current solution feasible

    int nRows() const { return rows.size(); }
    int nCols() const { return cols.size(); }

    static double DEtoBasisProb(double DE) { return DE<0.0?exp(DE):1.0; } // convert change in energy of a column to probability of proposal
//    static double DEtoBasisProb(double DE) { return exp(0.5*DE); } // convert change in energy of a column to probability of proposal

protected:
//    Proposal makeProposal();
    void insert(int i,int j,T entry) {
        cols[j].insert(i, entry);
        rows[i].insert(j, entry);
    }

    T infeasibility(int i, T Xi); // infeasibility factor for row i at value v

    double energy(int i, T Xi);

    void applyProposal(const Proposal &proposal, bool updateX);

    const std::vector<T> &nextSampleWithInfeasibleImportance();

//    const std::vector<T> &nextSample();
};


// copies the entries of 'tableau' into this, removing equality constraints
// and columns corresponding to basic vars
template<class T>
SparseBasisSampler<T>::SparseBasisSampler(
        const TableauNormMinimiser<T> &tableau,
        const std::vector<std::function<double(T)>> &tableaufactors,
        std::unique_ptr<PerturbableFunction<T,double>> &&importance,
        const std::vector<T> &initialState,
        double Kappa):
        kappa(Kappa),
        cols(tableau.nNonBasic()),
        rows(tableau.nAuxiliaryVars),
        basisDistribution(cols.size()),
        importanceFunc(std::move(importance)),
        currentDE(cols.size(),0.0),
        currentInfeasibility(0)
//        currentLogPiota(0.0)
{
    std::vector<int> tableauRowToBasisRow; // map from tableau row index to this row index, -1 means equality constraint
    tableauRowToBasisRow.reserve(tableau.rows.size());
    L.reserve(rows.size());
    U.reserve(rows.size());
    factors.reserve(rows.size());
    X.reserve(rows.size());
    delta.reserve(cols.size());
    for(int basisi=0, i=0; i<tableau.rows.size(); ++i) {
        if(tableau.isEqualityConstraint(i)) {
            tableauRowToBasisRow[i] = -1;
        } else {
            tableauRowToBasisRow[i] = basisi;
            rows[basisi].reserve(tableau.rows[i].size());
            L.push_back(tableau.L[i]);
            U.push_back(tableau.U[i]);
            factors.push_back(tableaufactors[i]);
            X.push_back(tableau.F[i]); // start with all X_N = 0
            ++basisi;
        }
    }
//    std::cout << "initial currentInfeasibility = " << currentInfeasibility << std::endl;
    assert(factors.size() == tableau.nAuxiliaryVars);

    // transfer tableau over to this
    for(int basisj=0, j=0; j<tableau.cols.size(); ++j) {
        if(!tableau.cols[j].isBasic) {
            cols[basisj].reserve(tableau.cols[j].size());
            delta.push_back((initialState.size()>j && initialState[j]>0)?-1:1);
            for(auto i : tableau.cols[j]) {
                int basisi = tableauRowToBasisRow[i];
                assert(basisi != -1);
                T Mij = tableau.rows[i].at(j);
                insert(basisi, basisj, Mij);
            }
            basisDistribution.set(basisj, DEtoBasisProb(currentDE[basisj]));
            ++basisj;
        }
    }
//    std::cout << "Initial basisDistribution" << basisDistribution << std::endl;
    assert(delta.size() == tableau.nNonBasic());

    // initialise X, currentE, currentInfeasibility and currentLogPiota
    currentE.reserve(rows.size());
    for(int basisi=0; basisi < rows.size(); ++basisi) {
        const SparseVec<T> &row = rows[basisi];
        for(int nzj=0; nzj < row.sparseSize(); ++nzj) {
            if(delta[row.indices[nzj]] == -1) X[basisi] += row.values[nzj];
        }
        currentE.push_back(energy(basisi, X[basisi]));
        currentInfeasibility += infeasibility(basisi, X[basisi]);
//        currentLogPiota += currentE[basisi];
//        std::cout << "Initial E[" << basisi << "] = " << currentE[basisi] << " " << X[basisi] << std::endl;
    }

    // initialise DE
    for(int basisj=0; basisj< cols.size(); ++basisj) {
        const SparseVec<T> &col = cols[basisj];
        for(int nzi = 0; nzi < col.sparseSize(); ++nzi) {
            int basisi = col.indices[nzi];
            currentDE[basisj] += energy(basisi, X[basisi] + delta[basisj] * col.values[nzi]) - currentE[basisi];
        }
//        std::cout << "Initial DE[" << basisj << "] = " << currentDE[basisj] << std::endl;
    }

    findInitialFeasibleSolution();

    // initialise importanceFunc and currentImportance
    importanceFunc->setState(X);
    currentLogImportance = importanceFunc->getLogValue(X);
}


template<class T>
SparseBasisSampler<T>::SparseBasisSampler(const WeightedFactoredConvexDistribution<T> &distribution, const std::vector<T> &initialState, double Kappa):
        SparseBasisSampler(
                TableauNormMinimiser<T>(distribution.constraints),
                distribution.factors,
                distribution.perturbableFunctionFactory(),
                initialState,
                Kappa)
{ }

template<class T>
SparseBasisSampler<T>::SparseBasisSampler(const WeightedFactoredConvexDistribution<T> &distribution, double Kappa):
        SparseBasisSampler(
                TableauNormMinimiser<T>(distribution.constraints),
                distribution.factors,
                distribution.perturbableFunctionFactory(),
                std::vector<T>(),
                Kappa)
{ }


//template<class T>
//void SparseBasisSampler<T>::initProbabilities() {
//    probContributions.resize(nCols());
//    for(int j=0; j<nCols(); ++j) {
//        probContributions[j].reserve(cols[j].sparseSize());
//        double Pj = 1.0;
//        for(int z=0; z < cols[j].sparseSize(); ++z) {
//            int i = cols[j].indices[z];
//            T Mij = cols[j].values[z];
//            double Xiprime = X[i] + Mij * Delta[j];
//            double contrib = std::min(1.0, marginalProbs(i, Xiprime)*iota(i, Xiprime)/(marginalProbs(i,X[i])*iota(i,X[i])));
//            probContributions[j].push_back(contrib);
//            Pj *= contrib;
//        }
//        basisDistribution.set(j, Pj);
//    }
//}


//template<class T>
//typename SparseBasisSampler<T>::Proposal SparseBasisSampler<T>::makeProposal() {
//    Proposal proposal(currentInfeasibility, currentLogPiota);
//
//    int j = basisDistribution(Random::gen); // choose basis to perturb
//    const SparseVec<T> &proposedBasis = cols[j];
//    proposal.changedIndices = proposedBasis.indices;
//    proposal.changedDE[j] = -currentDE[j];
//    T dj = Delta[j];
//    proposal.changedX.reserve(proposedBasis.sparseSize());
//    for(int nzi=0; nzi < proposedBasis.sparseSize(); ++nzi) {
//        int i = proposedBasis.indices[nzi]; // row that is affected by the change
//        T deltaXi = proposedBasis.values[nzi] * dj;
//        T newXi = X[i] + deltaXi; // the value of X[i] after swapping col j;
//        proposal.changedX.push_back(newXi);
//        double deltaEi = energy(i, newXi) - currentE[i]; // Energy of row i after flipping column j
//        proposal.logPiota += deltaEi;
//        proposal.infeasibility += infeasibility(i, newXi) - infeasibility(i, X[i]);
//        const SparseVec<T> &row = rows[i];
//        for (int nzj = 0; nzj < row.sparseSize(); ++nzj) {
//            int updatej = row.indices[nzj];
//            if(updatej != j) {
//                auto [changedDEIterator, wasInserted] = proposal.changedDE.try_emplace(updatej, currentDE[updatej]);
//                double deltaXiupdatej = row.values[nzj] * Delta[updatej]; // the change in X[i] on updating col updatej
//                changedDEIterator->second +=
//                        energy(i, newXi + deltaXiupdatej) - energy(i, X[i] + deltaXiupdatej) - deltaEi;
//            }
//        }
//    }
//
//    return proposal;
//}

// Take a sample
template<class T>
const std::vector<T> &SparseBasisSampler<T>::operator()() {
    return nextSampleWithInfeasibleImportance();
}

// Take a sample
//template<class T>
//const std::vector<T> &SparseBasisSampler<T>::nextSample() {
//
//    do {
////        std::cout << "Drawing proposal from distribution " << basisDistribution << "  " << currentDE << std::endl;
//
//        Proposal proposal(*this);
//
////        std::cout << "Got proposal col=" << proposal.proposedj << " inf=" << proposal.infeasibility << " DE=" << currentDE[proposal.proposedj] << std::endl;
//
//        if (proposal.infeasibility == 0) {       // --- transition to feasible state
//
//            importanceFunc->perturbWithUndo(X, proposal.changedXIndices);
//            for (int nzi = 0; nzi < proposal.changedXIndices.size(); ++nzi) {
//                // update X and transform changedX into undo
//                std::swap(X[proposal.changedXIndices[nzi]], proposal.changedXValues[nzi]);
//            }
//
//            double proposedLogImportance = importanceFunc->getLogValue(X);
////            std::cout << "Proposed importance = " << proposedImportance << std::endl;
//            double ratioOfSums = proposal.calcRatioOfSums(basisDistribution);
//            double acceptance = exp(proposedLogImportance - currentLogImportance)*ratioOfSums;
////            std::cout << "Ratio of sums = " << ratioOfSums << " acceptance = " << acceptance << std::endl;
//
//            if (Random::nextDouble() < acceptance) {     // --- accept ---
////                std::cout << "Accepting feasible" << std::endl;
//                stats.addSample(true, currentInfeasibility == 0, true);
//                currentLogImportance = proposedLogImportance;
//                applyProposal(proposal, false);
//            } else {                                    // --- reject ---
////                std::cout << "rejecting feasible" << std::endl;
//                // undo changes to X
//                for (int nzi = 0; nzi < proposal.changedXIndices.size(); ++nzi) {
//                    X[proposal.changedXIndices[nzi]] = proposal.changedXValues[nzi];
//                }
//                // undo changes to importanceFunc
//                importanceFunc->undoLastPerturbation();
//                stats.addSample(false, currentInfeasibility == 0, true);
//            }
//        } else {                                // --- transition to infeasible state
//            double logRatioOfSums = proposal.calcLogRatioOfSums(basisDistribution);
//            double logAcceptance = logRatioOfSums - currentLogImportance;
////            std::cout << "Ratio of sums = " << ratioOfSums << " acceptance = " << acceptance << std::endl;
//            if (Random::nextDouble() < exp(logAcceptance)) {     // --- accept ---
////                std::cout << "Accepting infeasible proposal inf=" << proposal.infeasibility << std::endl;
//                stats.addSample(true, currentInfeasibility == 0, false);
//                currentLogImportance = 0.0;
//                applyProposal(proposal, true);
//                importanceFunc->perturb(X, proposal.changedXIndices);
//            } else {                                    // --- reject ---
////                std::cout << "rejecting infeasible proposal inf=" << proposal.infeasibility << std::endl;
//                stats.addSample(false, currentInfeasibility == 0, false);
//            }
//        }
//    } while(currentInfeasibility > 0);
//    return X;
//}


// Calculate importance for infeasible states as well as feasible
template<class T>
const std::vector<T> &SparseBasisSampler<T>::nextSampleWithInfeasibleImportance() {

    do {
//        std::cout << "Drawing proposal from distribution " << basisDistribution << "  " << currentDE << std::endl;
        Proposal proposal(*this);

//        std::cout << "Got proposal col=" << proposal.proposedj << " inf=" << proposal.infeasibility << " DE=" << currentDE[proposal.proposedj] << std::endl;

        importanceFunc->perturb(X, proposal.changedXIndices);
        for (int nzi = 0; nzi < proposal.changedXIndices.size(); ++nzi) {
            // update X and transform changedX into undo
            std::swap(X[proposal.changedXIndices[nzi]], proposal.changedXValues[nzi]);
        }

        double proposedLogImportance = importanceFunc->getLogValue(X);
//            std::cout << "Proposed importance = " << proposedImportance << std::endl;
        double ratioOfSums = proposal.calcRatioOfSums(basisDistribution);
        double acceptance = exp(proposedLogImportance - currentLogImportance)*ratioOfSums;

//            std::cout << "Ratio of sums = " << ratioOfSums << " acceptance = " << acceptance << std::endl;

        if (Random::nextDouble() < acceptance) {     // --- accept ---
//                std::cout << "Accepting feasible" << std::endl;
            stats.addSample(true, currentInfeasibility == 0, proposal.infeasibility == 0);
            currentLogImportance = proposedLogImportance;
            applyProposal(proposal, false);
        } else {                                    // --- reject ---
//                std::cout << "rejecting feasible" << std::endl;
            // undo changes to X
            for (int nzi = 0; nzi < proposal.changedXIndices.size(); ++nzi) {
                X[proposal.changedXIndices[nzi]] = proposal.changedXValues[nzi];
            }
            // undo changes to importanceFunc
            importanceFunc->undo();
            stats.addSample(false, currentInfeasibility == 0, proposal.infeasibility == 0);
        }
    } while (currentInfeasibility > 0);
    return X;
}


// Makes the current solution feasible by choosing pivots from the raw column probability
template<class T>
void SparseBasisSampler<T>::findInitialFeasibleSolution() {
    int iterations = 0;

//    for(int burnin=0; burnin < 1000000; ++burnin) {  // TODO: JUST A TEST!!!!!!!!!!!!!!!
//        Proposal proposal(*this);
//        applyProposal(proposal, true);
//    }

    while(currentInfeasibility > 0) {
        Proposal proposal(*this);
        applyProposal(proposal, true);
        ++iterations;
    }
    std::cout << "Found feasible solution in " << iterations << " iterations" << std::endl;
}


template<class T>
T SparseBasisSampler<T>::infeasibility(int i, T Xi) {
    if(Xi < L[i]) return L[i]-Xi;
    if(Xi > U[i]) return Xi-U[i];
    return 0;
}



template<class T>
double SparseBasisSampler<T>::Proposal::calcRatioOfSums(const MutableCategoricalArray &basisDistribution) {
    double changeInSum = 0.0;
    for(const auto [changedj, newDeltaEj]: changedDE) {
//        std::cout << "newDeltaEj = " << newDeltaEj << std::endl;
        double oldPj = basisDistribution[changedj];
        changeInSum += SparseBasisSampler<T>::DEtoBasisProb(newDeltaEj) - oldPj;
    }
//    std::cout << "basis sum of probs = " << basisDistribution.sum() << std::endl;
    double ratio = basisDistribution.sum()/(basisDistribution.sum() + changeInSum);
//    std::cout << ratio << std::endl;
    return ratio;
}


template<class T>
SparseBasisSampler<T>::Proposal::Proposal(SparseBasisSampler<T> &parent):
proposedj(parent.basisDistribution(Random::gen)),
infeasibility(parent.currentInfeasibility),
//logPiota(parent.currentLogPiota),
changedXIndices(parent.cols[proposedj].indices)
{
    const SparseVec<T> &proposedBasis = parent.cols[proposedj];
    changedDE[proposedj] = -parent.currentDE[proposedj];
    T dj = parent.delta[proposedj];
    changedXValues.reserve(proposedBasis.sparseSize());
    for(int nzi=0; nzi < proposedBasis.sparseSize(); ++nzi) {
        int i = proposedBasis.indices[nzi]; // row that is affected by the change
        T deltaXi = proposedBasis.values[nzi] * dj;
        T newXi = parent.X[i] + deltaXi; // the value of X[i] after swapping col j;
        changedXValues.push_back(newXi);
        double deltaEi = parent.energy(i, newXi) - parent.currentE[i]; // change in energy of row i by flipping column j
//        logPiota += deltaEi;
        infeasibility += parent.infeasibility(i, newXi) - parent.infeasibility(i, parent.X[i]);
        const SparseVec<T> &row = parent.rows[i];
        for (int nzj = 0; nzj < row.sparseSize(); ++nzj) {
            int updatej = row.indices[nzj];
            if(updatej != proposedj) {
                auto [changedDEIterator, wasInserted] = changedDE.try_emplace(updatej, parent.currentDE[updatej]);
                double deltaXiupdatej = row.values[nzj] * parent.delta[updatej]; // the change in X[i] on updating col updatej
                changedDEIterator->second +=
                        parent.energy(i, newXi + deltaXiupdatej) - parent.energy(i, parent.X[i] + deltaXiupdatej) - deltaEi;
            }
        }
    }
}


// assuming the last proposal was accepted, beta has been
// updated accordingly and changedDeltaE constains the changes
// to deltaE caused by the pivot.
template<class T>
void SparseBasisSampler<T>::applyProposal(const Proposal &proposal, bool updateX) {
//    currentLogPiota = proposal.logPiota;
    currentInfeasibility = proposal.infeasibility;

    for(int nzi=0; nzi<proposal.changedXIndices.size(); ++nzi) {
        int i = proposal.changedXIndices[nzi];
        if(updateX) X[i] = proposal.changedXValues[nzi];
        currentE[i] = energy(i, X[i]);
    }

    for(const auto [changedj, newDeltaEj]: proposal.changedDE) {
        currentDE[changedj] = newDeltaEj;
        basisDistribution[changedj] = DEtoBasisProb(newDeltaEj);
    }
    delta[proposal.proposedj] *= -1;
}


// energy of the i'th row, for a given X[i]
template<class T>
double SparseBasisSampler<T>::energy(int i, T Xi) {
    if(Xi < L[i]) return factors[i](L[i]) - kappa*(L[i] - Xi);
    if(Xi > U[i]) return  factors[i](U[i]) - kappa*(Xi - U[i]);
    return factors[i](Xi);
}


template<class T>
std::ostream &operator <<(std::ostream &out, const SparseBasisSampler<T> &basis) {
    for(int i=0; i<basis.nRows(); ++i) {
        out << basis.L[i] << "\t<=\t" << basis.X[i] << "\t=\t";
        std::vector<T> row = basis.rows[i].toDense(basis.nCols());
        for(int j=0; j<basis.nCols(); ++j) {
            out << row[j] << "\t";
        }
        out << "<=\t" << basis.U[i] << "\t["
            << exp(basis.factors[i](basis.L[i])) << " ... "
            << exp(basis.factors[i](basis.U[i])) << "]" << std::endl;
    }
    return out;
}



#endif //ABMCMC_SPARSEBASISSAMPLER_H
