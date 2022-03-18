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

#include <vector>
#include "TableauNormMinimiser.h"
#include "MutableCategoricalArray.h"


template<class T>
class SparseBasisSampler {
public:
    // tableau
    std::vector<SparseVec<T>> cols; // the basis vectors
    std::vector<SparseVec<T>> rows; //
    std::vector<T>      H; // The upper bounds of the hyperrectangle by row (see intro notes above)
    std::vector<double>   marginalGradE; // Rate of change of energy with X[i], by i.
    double              kappa; // rate of decay of probability outside the hyper-rectangle

    std::vector<T>      X; // the current point on the lattice by row.
    std::vector<T>      Delta; // direction of perturbation of each basis by column (assume act Fermionicity)
    MutableCategoricalArray basisDistribution; // distribution of proposing to update basis by column
    std::vector<double> currentDE; // change in energy by flipping the value of the j'th column
    std::vector<double> currentE;  // current energy of i'th row
    T                   currentInfeasibility;
    double              currentLogPiota; // log of marginalised probability
    double              currentImportance; // real prob over estimated prob of current state (=1 if infeasible)

    class Proposal {
    public:
        Proposal(T initialInfeasibility, double initialLogPiota): infeasibility(initialInfeasibility), logPiota(initialLogPiota) {}

        int                     j;              // col to update
        SparseVec<T>            changedX;       // changed elements of X
        std::map<int,double>    changedDE;      // change in Delta-Energy by column
        T                       infeasibility;
        double                  logPiota;

        double calcRatioOfSums(const MutableCategoricalArray &basisDistribution);
    };

    SparseBasisSampler(TableauNormMinimiser &tableau, const std::vector<double> & marginalProbs);

    const std::vector<T> nextSample();

    int nRows() const { return H.size(); }
    int nCols() const { return this->size(); }

protected:
    Proposal makeProposal();
    void insert(int i,int j,T entry) {
        cols[j].insert(i, entry);
        rows[i].insert(j, entry);
    }

    double infeasibility(int i, T v); // infeasibility factor for row i at value v

    double energy(int i, double b);

    void applyProposal(const Proposal &proposal);

};


// marginalProbs takes an index into X and a value for X_i and returns the marginal prob of X_i having that value.
template<class T>
SparseBasisSampler<T>::SparseBasisSampler(TableauNormMinimiser &tableau, const std::vector<double> &marginalGradientE):
        cols(tableau.nNonBasic()),
        rows(tableau.cols.size()+tableau.nAuxiliaryVars),
        H(rows.size()),
        X(rows.size(), 0),
        Delta(tableau.nNonBasic(), 1),
        basisDistribution(tableau.nNonBasic()),
        marginalGradE(marginalGradientE)
{
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

    for(int i=0; i < tableau.Hc.size(); ++i) H[i] = tableau.Hc[i];
    for(int a=0; a < tableau.nAuxiliaryVars; ++a) H[a+tableau.Hc.size()] = tableau.Ha[a+1];

    for(int i=0; i < tableau.basis.size(); ++i) {
        int basisi = tableau.basis[i]<0?tableau.cols.size()-1-tableau.basis[i]:tableau.basis[i];
        X[basisi] = tableau.F[i];
    }
//    initProbabilities();
//    initIntersections();
}

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


template<class T>
typename SparseBasisSampler<T>::Proposal SparseBasisSampler<T>::makeProposal() {
    Proposal proposal;

    proposal.j = basisDistribution(Random::gen);
    const SparseVec<T> &col = cols[proposal.j];
    proposal.changedDE[proposal.j] = -currentDE[proposal.j];
    T dj = Delta[proposal.j];
    for(int nzi=0; nzi<col.sparseSize(); ++nzi) {
        int i = col.indices[nzi]; // row that is affected by the change
        T deltaXi = col.values[nzi] * dj;
        T newXi = X[i] + deltaXi; // the value of X[i] after swapping col j;
        proposal.changedX.insert(i, newXi);
        double deltaEi = energy(i, newXi) - currentE[i]; // Energy of row i after flipping column j
        proposal.logPiota += deltaEi;
        proposal.infeasibility += infeasibility(i, newXi) - infeasibility(i, X[i]);
        const SparseVec<T> &row = rows[i];
        for (int nzj = 0; nzj < row.sparseSize(); ++nzj) {
            int updatej = row.indices[nzj];
            if(updatej != proposal.j) {
                auto [proposedDEelement, wasInserted] = proposal.changedDE.try_emplace(updatej, currentDE[updatej]);
                double deltaXiupdatej = row.values[nzj] * Delta[updatej]; // the change in X[i] on updating col updatej
                proposedDEelement->second +=
                        energy(i, newXi + deltaXiupdatej) - energy(i, X[i] + deltaXiupdatej) - deltaEi;
            }
        }
    }

    return proposal;
}


template<class T>
const std::vector<T> SparseBasisSampler<T>::nextSample() {

    Proposal proposal = makeProposal();

//    const SparseVec<T> &basis = cols[j];
//    // update X
//    SparseVec<T> deltaX;
//    for(int z=0; z<basis.sparseSize(); ++z) {
//        deltaX.insert(basis.indices[z], Delta[j]*basis.values[z]);
////        X[basis.indices[z]] += Delta[j]*basis.values[z];
//    }
//
//    std::map<int,double> changedDeltaE = calcDeltaEChanges(j);
//    double ratioOfSums = calcRatioOfSums(j, changedDeltaE);

    double proposedImportance = (proposal.infeasibility == 0?
                    exp(logP(proposal.changedX) - proposal.logPiota)
                    :1.0
            );
    double acceptance = proposal.calcRatioOfSums(basisDistribution)*currentImportance/proposedImportance;

    if(Random::nextDouble() < acceptance) {
        // accept
        currentImportance = proposedImportance;
        applyProposal(proposal);
    } else {
        // reject
    }

    std::vector<int> touchedCols; // can be pre-computed for each j
    // TODO: for each touched factor, we need to know column (j) and sparse entry (z),
    //  so could just store a vector of touched <j,z> pairs for each column
    //  or even better a vector of <j,<z_1..z_n>> pairs where we aggregate over all paris with the same j
    return X;
}

template<class T>
double SparseBasisSampler<T>::infeasibility(int i, T v) {
    if(v < 0) return -v;
    if(v > H[i]) return (H[i]-v);
    return 0;
}



template<class T>
double SparseBasisSampler<T>::Proposal::calcRatioOfSums(const MutableCategoricalArray &basisDistribution) {
    double changeInSum = 0.0;
    for(const auto [changedj, newDeltaEj]: changedDE) {
        double oldPj = basisDistribution[changedj];
        changeInSum += (newDeltaEj>0.0?exp(newDeltaEj):1.0) - oldPj;
    }
    double ratio = basisDistribution.sum()/(basisDistribution.sum() + changeInSum);
//    std::cout << ratio << std::endl;
    return ratio;
}



// assuming the last proposal was accepted, beta has been
// updated accordingly and changedDeltaE constains the changes
// to deltaE caused by the pivot.
template<class T>
void SparseBasisSampler<T>::applyProposal(const Proposal &proposal) {
//    const SparseVec<T> &col = cols[proposal.j];
//    for(int nzi=0; nzi<col.sparseSize(); ++nzi) {
//        int i = col.indices[nzi];
//        currentE[i] = energy(i, X[i]);
//    }
    currentLogPiota = proposal.logPiota;
    currentInfeasibility = proposal.infeasibility;
    for(int nzi=0; nzi<proposal.changedX.sparseSize(); ++nzi) {
        int i = proposal.changedX.indices[nzi];
        X[i] = proposal.changedX.values[nzi];
        currentE[i] = energy(i, X[i]);
    }

    for(const auto [changedj, newDeltaEj]: proposal.changedDE) {
        currentDE[changedj] = newDeltaEj;
        basisDistribution[changedj] = newDeltaEj>0.0?exp(newDeltaEj):1.0;
    }
}


// energy of the i'th row, for a given X[i]
template<class T>
double SparseBasisSampler<T>::energy(int i, double Xi) {
    if(Xi > H[i]) return kappa*(Xi - H[i]) + marginalGradE[i]*H[i];
    if(Xi < 0) return -kappa*Xi;
    return Xi * marginalGradE[i];
}


template<class T>
std::ostream &operator <<(std::ostream &out, const SparseBasisSampler<T> &basis) {
    for(int i=0; i<basis.nRows(); ++i) {
        out << basis.X[i] << " = ";
        for(int j=0; j<basis.nCols(); ++j) {
            out << basis[j][i] << "\t";
        }
        out << " <= " << basis.H[i] << std::endl;
    }
    return out;
}



#endif //ABMCMC_SPARSEBASISSAMPLER_H
