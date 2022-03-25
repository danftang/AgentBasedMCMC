// Samples from a lattice of points in a plane represented as a set of basis vectors, stored in sparse format.
// The valid points of the lattice are contained within a hyperrectangle with one corner at the origin
// and a diagonally opposite corner at "H", so that the valid points can be defined as the set
//
// { X : X = BX' + F, 0 <= X <= H }
//
// where the columns of B are the basis vectors, and the elements of $X'$ are integers.
//
// The probability of each point must be expressible as a product of unnormalised marginals
// P(X) = A\prod_i P_i(X_i) \prod_i I(X_i)
// where (since we're dealing only with the act Fermionic case) we assume P_i takes the form
//
// P_i(X_i) = 1 if X_i <= 0 or pi_i if X >= 1
// for some constant pi_i. So P_i(X'_i)/P_i(X_i) is always either pi_i or 1/pi_i.
//
// and I(X_i) = 1 if 0 <= X_i <= H_i, e^{k X_i} if X_i < 0 or e^{k(H_i-Xi)} if X_i > H_i
//
// Created by daniel on 07/03/2022.
//

#ifndef ABMCMC_SPARSEBASISSAMPLER_H
#define ABMCMC_SPARSEBASISSAMPLER_H

#include <utility>
#include <vector>
#include "TableauNormMinimiser.h"
#include "MutableCategoricalArray.h"

// TODO: need to generalise marginalGradE to general function, get rid of X (now in importance function),
//  integrate importance function and finish construction from WeightedFactoredConvexDistribution.

template<class T>
class SparseBasisSampler {
public:
    static constexpr double   kappa = 8.0; // rate of decay of probability outside the hyper-rectangle

    // tableau
    std::vector<SparseVec<T>>   cols;   // the basis vectors
    std::vector<SparseVec<T>>   rows;   //
    std::vector<T>              H;      // The upper bounds of the hyperrectangle by row (see intro notes above)
    std::vector<std::function<double(T)>> factors;

    std::vector<T>              X;      // the current point on the lattice by row.
    std::vector<T>              delta;  // direction of perturbation of each basis by column (assume act Fermionicity)
    MutableCategoricalArray     basisDistribution; // distribution of proposing to update basis by column
    std::vector<double>         currentDE;      // change in energy by flipping the value of the j'th column
    std::vector<double>         currentE;       // current energy of i'th row
    T                           currentInfeasibility;
    double                      currentLogPiota; // log of marginalised probability, equal to minus the total current energy
    double                      currentImportance; // real prob over estimated prob of current state (=1 if infeasible)

    std::unique_ptr<PerturbableFunction<T,double>>    importanceFunc;

    class Proposal {
    public:
        Proposal(SparseBasisSampler<T> &parent);

        int proposedj;                           // column index of proposed basis
        const std::vector<int> &changedXIndices; // indices of X that are changed in the proposal
        std::vector<T>          changedXValues;  // new values of elements of X, by sparse index of j'th basis
        std::map<int,double>    changedDE;       // change in Delta-Energy by column
        T                       infeasibility;   // proposed infeasibility
        double                  logPiota;        // proposed log P_i

        double calcRatioOfSums(const MutableCategoricalArray &basisDistribution);
    };

    SparseBasisSampler(const TableauNormMinimiser<T> &tableau, const std::vector<std::function<double(T)>> &tableauFactors, std::unique_ptr<PerturbableFunction<T,double>> &&importanceFunc);

    explicit SparseBasisSampler(const WeightedFactoredConvexDistribution<T> &distribution);

    const std::vector<T> nextSample();

    int nRows() const { return rows.size(); }
    int nCols() const { return cols.size(); }

protected:
//    Proposal makeProposal();
    void insert(int i,int j,T entry) {
        cols[j].insert(i, entry);
        rows[i].insert(j, entry);
    }

    T infeasibility(int i, T Xi); // infeasibility factor for row i at value v

    double energy(int i, T Xi);

    void applyProposal(const Proposal &proposal, bool updateX);

};


// marginalProbs takes an index into X and a value for X_i and returns the marginal prob of X_i having that value.
template<class T>
SparseBasisSampler<T>::SparseBasisSampler(
        const TableauNormMinimiser<T> &tableau,
        const std::vector<std::function<double(T)>> &tableaufactors,
        std::unique_ptr<PerturbableFunction<T,double>> &&importance):
        cols(tableau.nNonBasic()),
        rows(tableau.cols.size()+tableau.nAuxiliaryVars),
        H(rows.size()),
        X(rows.size(), 0),
        delta(cols.size(), 1),
        basisDistribution(cols.size()),
        factors(rows.size()),
        importanceFunc(std::move(importance)),
        currentDE(cols.size()),
        currentE(X.size())
{
    // transfer tableau over to this
    int basisj = 0;
    for(int j=0; j<tableau.cols.size(); ++j) {
        if(!tableau.cols[j].isBasic) {
            cols[basisj].reserve(tableau.cols[j].size() + 1);
            for(auto i : tableau.cols[j]) {
                int basisi = tableau.basis[i];
                if(basisi < 0) basisi = -basisi + tableau.cols.size() - 1; // transform auxiliaries to after end of X
                insert(basisi, basisj, tableau.rows[i].at(j));
            }
            insert(j, basisj, 1);
            ++basisj;
        }
    }

    // initialise H and factors
    for(int j=0; j < tableau.cols.size(); ++j) {
        H[j] = tableau.Hc[j];
        factors[j] = tableaufactors[tableau.constraintIndexByCol[j]];
        std::cout << "Setting up factor " << j << " -> " << factors[j](0) << ", " << factors[j](1) << std::endl;
    }
    for(int a=0; a < tableau.nAuxiliaryVars; ++a) {
        int j = a+tableau.cols.size();
        H[j] = tableau.Ha[a+1];
        factors[j] = tableaufactors[tableau.constraintIndexByAuxiliary[a+1]];
        std::cout << "Setting up factor " << j << " -> " << factors[j](0) << ", " << factors[j](1) << std::endl;
    }

    // initialise X
    for(int i=0; i < tableau.basis.size(); ++i) {
        int basisi = tableau.basis[i]<0?tableau.cols.size()-1-tableau.basis[i]:tableau.basis[i];
        X[basisi] = tableau.F[i];
    }





    // initialise currentE and currentInfeasibility
    currentInfeasibility = 0;
    currentLogPiota = 0.0;
    for(int i=0; i<rows.size(); ++i) {
        currentE[i] = energy(i, X[i]);
        currentInfeasibility += infeasibility(i, X[i]);
        currentLogPiota -= currentE[i];
        std::cout << "Initial E[" << i << "] = " << currentE[i] << " " << X[i] << std::endl;
    }

    // initialise currentDE and basisDistribution
    for(int j=0; j<cols.size(); ++j) {
        currentDE[j] = 0.0;
        for(int nzi = 0; nzi < cols[j].sparseSize(); ++nzi) {
            int changedRow = cols[j].indices[nzi];
            currentDE[j] += energy(changedRow, X[changedRow] + delta[changedRow] * cols[j].values[nzi]) - currentE[changedRow];
        }
        basisDistribution.set(j, currentDE[j]>0.0?exp(-currentDE[j]):1.0);
        std::cout << "Initial DE[" << j << "] = " << currentDE[j] << std::endl;
    }

    // initialise importanceFunc and currentImportance
    importanceFunc->setState(X);
    if(currentInfeasibility == 0) {
        currentImportance = importanceFunc->getValue(X);
    } else {
        currentImportance = 1.0;
    }

}


template<class T>
SparseBasisSampler<T>::SparseBasisSampler(const WeightedFactoredConvexDistribution<T> &distribution):
        SparseBasisSampler(TableauNormMinimiser<T>(distribution.constraints), distribution.factors, distribution.perturbableFunctionFactory())
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


template<class T>
const std::vector<T> SparseBasisSampler<T>::nextSample() {

    Proposal proposal(*this);

    std::cout << "Got proposal col=" << proposal.proposedj << " inf=" << proposal.infeasibility << " DE=" << currentDE[proposal.proposedj] << std::endl;

    if(proposal.infeasibility == 0) {       // --- transition to feasible state

        importanceFunc->perturbWithUndo(X, proposal.changedXIndices);
        for(int nzi=0; nzi < proposal.changedXIndices.size(); ++nzi) {
            // update X and transform changedX into undo
            std::swap(X[proposal.changedXIndices[nzi]], proposal.changedXValues[nzi]);
        }

        double proposedImportance = importanceFunc->getValue(X);
        double ratioOfSums = proposal.calcRatioOfSums(basisDistribution);
        double acceptance = ratioOfSums*currentImportance/proposedImportance;

        std::cout << "Ratio of sums = " << ratioOfSums << " acceptance = " << acceptance << std::endl;

        if(Random::nextDouble() < acceptance) {     // --- accept ---
            currentImportance = proposedImportance;
            applyProposal(proposal, false);
        } else {                                    // --- reject ---
            std::cout << "rejecting feasible" << std::endl;
            // undo changes to X
            for(int nzi=0; nzi < proposal.changedXIndices.size(); ++nzi) {
                X[proposal.changedXIndices[nzi]] = proposal.changedXValues[nzi];
            }
            // undo changes to importanceFunc
            importanceFunc->undoLastPerturbation();
        }
    } else {                                // --- transition to infeasible state
        double acceptance = proposal.calcRatioOfSums(basisDistribution)*currentImportance;
        if(Random::nextDouble() < acceptance) {     // --- accept ---
            currentImportance = 1.0;
            applyProposal(proposal, true);
            importanceFunc->perturb(X, proposal.changedXIndices);
        } else {
            std::cout << "rejecting infeasible" << std::endl;
        }
    }

    return X;
}

template<class T>
T SparseBasisSampler<T>::infeasibility(int i, T Xi) {
    if(Xi < 0) return -Xi;
    if(Xi > H[i]) return Xi-H[i];
    return 0;
}



template<class T>
double SparseBasisSampler<T>::Proposal::calcRatioOfSums(const MutableCategoricalArray &basisDistribution) {
    double changeInSum = 0.0;
    for(const auto [changedj, newDeltaEj]: changedDE) {
        double oldPj = basisDistribution[changedj];
        changeInSum += (newDeltaEj>0.0?exp(-newDeltaEj):1.0) - oldPj;
    }
    double ratio = basisDistribution.sum()/(basisDistribution.sum() + changeInSum);
//    std::cout << ratio << std::endl;
    return ratio;
}


template<class T>
SparseBasisSampler<T>::Proposal::Proposal(SparseBasisSampler<T> &parent):
proposedj(parent.basisDistribution(Random::gen)),
infeasibility(parent.currentInfeasibility),
logPiota(parent.currentLogPiota),
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
        logPiota += deltaEi;
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
    currentLogPiota = proposal.logPiota;
    currentInfeasibility = proposal.infeasibility;

    for(int nzi=0; nzi<proposal.changedXIndices.size(); ++nzi) {
        int i = proposal.changedXIndices[nzi];
        if(updateX) X[i] = proposal.changedXValues[nzi];
        currentE[i] = energy(i, X[i]);
    }

    for(const auto [changedj, newDeltaEj]: proposal.changedDE) {
        currentDE[changedj] = newDeltaEj;
        basisDistribution[changedj] = newDeltaEj>0.0?exp(-newDeltaEj):1.0;
    }
    delta[proposal.proposedj] *= -1;
}


// energy of the i'th row, for a given X[i]
template<class T>
double SparseBasisSampler<T>::energy(int i, T Xi) {
    if(Xi > H[i]) return  factors[i](H[i]) + kappa*(Xi - H[i]);
    if(Xi < 0) return factors[i](0) - kappa*Xi;
    return factors[i](Xi);
}


template<class T>
std::ostream &operator <<(std::ostream &out, const SparseBasisSampler<T> &basis) {
    for(int i=0; i<basis.nRows(); ++i) {
        out << basis.X[i] << " = ";
        std::vector<T> row = basis.rows[i].toDense(basis.nCols());
        for(int j=0; j<basis.nCols(); ++j) {
            out << row[j] << "\t";
        }
        out << " <= " << basis.H[i] << std::endl;
    }
    return out;
}



#endif //ABMCMC_SPARSEBASISSAMPLER_H
