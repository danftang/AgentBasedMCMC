
// Find basis that approximately maximises the Kolmogorov-Smirnov entropy
// of the resulting Markov process.
//
// Kolmogorov-Smirnov entropy of a Markov process is
//
// S = -\sum_x P(x) \sum_y P(x->y)log(P(x->y))
//
// where P(x) is the stationary probability of being in state x
// and P(x->y) is the probability of transitioning from x to y.
// Here we don't include loop states in the sum so as to penalise against
// rejection.
//
// Notice that this is the expectation over the stationary distribution of the entropy of the transition distribution
//
// Given that the tableau will be used to create a set of basis vectors, b_j, from the columns to
// make a Markov process
//
// P(x->x+b_j) = (P(x+b_j)/P(x))^0.5 / s_x
// and
// P(x+b_j->x) = (P(x)/P(x+b_j))^0.5 / s_{x+b_j}
// where
// s_x = (\sum_z (P(z)/P(x))^0.5)
//
// then, considering pairs of transitions x->y and y->x
//
// 2S = -\sum_x  \sum_j P(x)((P(x+b_j)/P(x))^0.5 / s_x) log((P(x+b_j)/P(x))^0.5 / s_x) + P(x+b_j)((P(x)/P(x+b_j))^0.5 / s_{x+b_j}) log((P(x)/P(x+b_j))^0.5 / s_y)
//    = \sum_{x,j} (P(x)P(x+b_j))^0.5 * (0.5log(P(x+b_j)/P(x))/s_x + 0.5log(P(x)/P(x+b_j))/s_{x+b_j} - log(s_x)/s_x - log(s_{x+b_j})/s_{x+b_j})
// where
// If we make the approximation that s_x = s_{x+b_j} = s, this simplifies to
//
// S = (log(s)/s) * \sum_{x,j} (P(x)P(x+b_j))^0.5
// where
// b_j is the j'th basis
//
// But if we approximate s = \sum_{x,j} (P(x)P(x+b_j))^0.5 (i.e. E_P(x)[s_x]) then
//
// S = log(s) so maximising s maximises S.
//
// If we assume that P(x) = \prod_k f_k(x_k), where x_k is a single element if x then,
//
// s = \sum_j S_j
// where
// S_j = \prod_k F_k(b_{jk})            [notice that S_k(0) = 1]
//   is the expected entropy of the j'th column and
// F_k(d) = \sum_{x} (f_k(x)f_k(x+d))^0.5
//   is the expected entropy of a perturbation of d in the k'th factor
// b_{jk} is the k'th element of the j'th basis vector
// ...or if we want to get the expectation over a non-stationary prob dist use
// F_k(d) = \sum_{x} 0.5*P(x)(f_k(x+d)/f_k(x))^0.5 + 0.5*P(x+d)(f_k(x)/f_k(x+d))^0.5 ...?
//
// On pivot we set the row to zero.
// We approximate the change in entropy on pivoting at point (i,j) by assuming no overlap of elements so that
//
// DS_ij = (R'_i-S_j/F_i(M_ij))*S_j/F_i(M_ij) - R_i    [what if an element of row i is non-unit? Assume entropy contrib falls to zero, so contrib to R'_i on that col goes to zero?]
//
// where
// R_i = \sum_{j \in Row_i} S_j
//   is the sum of column entropies in row i
// R'_i = \sum_{j in Row_i} S_j/F_i(M_ij)
//   is the sum of column entropies after deletion of row i
//
// ...or...total entropy after a pivot at i,j:
//
// S_ij = \sum_{l in Row_i l != j} ( S'_ij S'_il - S_l )
//
// where S'_ij = S_j/F_i(M_ij)
//
// ...or... should we just calculate the absolute entropy after pivot? That way we only need to
// store 1 R''_i (S if row i was to have entropy of 1)
// entropy after pivot at i,j assuming no overlap
// s_ij = \sum_l S_l

// ...or... we could restrict ourselves to pivots that increase entropy, this only happens when
// we have a column with an element whose function probability is less than the product of
// all other elements in the column
// what's more, an entropy increasing pivot can only create more entropy increasing pivots (not
// invalidate existing ones) so it doesn't matter which order we do them in! [unless there's more
// than one entropy increasing pivot on a single row...so just choose the best per column then per row]
// if there's no entropy increasing pivots, then widen?

// F_i can be pre-calculated for 1, 2 and 3
// We store current values of S_j, R_i and R'_i, and the best pivot on each row
// On pivot at (i,j), we re-calculate S_k for each non-zero k in row i
// incrementing the suns R_i and R'_i for each entry in each updated column
//
// The best pivot column on a given row is the one that minimises
//
// |S_j/F_i(M_ij) - 0.5*R'_i|
//
// This is always the largest S_j/F_i(M_ij). So, we only need to update
// best pivot if a column that currently has a best pivot reduces in prob
// (need to re-check row for best pivot) or a column increases prob (need to
// compare each element in column with current best pivot of element row),
// or on rows where there is fill-in (need to check against current best pivot on
// filled-in row)
//
// Pivot(i,j):
// 0) For each non-zero column, l, on pivot row i,
//      for each non-zero row, k, on column l
//          subtract contribs from R_i and R'_i
//          add to set (vector) of affected rows (extract node?)
//          if this element is row-best-pivot, invalidate (and add to vector of invalidated rows?)
// 1) Do pivot on elements,
// 2) For each non-zero column, l, on pivot row i,
//      update S_l
//      for each non-zero row, k, on column l and not in column j
//          add contrib to R_k and R'_k
//          if row-best-pivot is valid, compare against current row-best-pivot and update if better
// 3) re-calculate invalid row-best-pivots
// 4) for each affected row
//      update priority queue entry of this row
//
//
//
// ... or ...
// We can consider the entropy in function space S(DF) = \sum_X (\prod_k F(X_k)F(X_k+DF_k))^0.5
// which has a maximum at DF = 0 (and rate of change 0) so we can approximate this
// by the curvature at 0. but the sum is over X in function space, not trajectory space
// In this space, there is no cross-correlation, so the curvature can be completely describeed
// by the curvature in each individual dimension (assuming monovariate functions)
//
// If we assume that the functions are all monovariate and normalised, then their product
// makes a distribution over "function space". If we augment the function space with a
// dimension for each constraint then we have a distribution over augmented function space
// and the function coeffs and constraints define a set of basis vectors in this augmented
// space. The sum over X is the sum over the augmented distribution, without any constraints.
//...or we can approximate the constrained distribution with any distribution that is fully
// factorised in the function space. Or, if we make auxiliary vars for each factor then the
// function space becomes the augmented trajectory space...
//
// if f(x) = e^k|x|
// \int e^k|x|e^k(|x-Dx|) dx = \int_-inf^0 e^-2kx e^kDx + \int_0^Dx e^kDx + \int_Dx^inf e^2kx e^-kDx
// = -e^kDx/2k + Dxe^Dx + -e^kDx/2k
// =  (Dx - 1/k)e^kDx
//
// Created by daniel on 17/10/22.

#ifndef ABMCMC_TABLEAUENTROPYMAXIMISER_H
#define ABMCMC_TABLEAUENTROPYMAXIMISER_H

#include <set>
#include <list>
#include <map>
#include "include/SparseVec.h"
#include "ConstrainedFactorisedDistribution.h"
#include "include/SortedSparseVec.h"
#include "include/StlStream.h"

template<class COEFF>
class TableauEntropyMaximiser {
public:
    static constexpr int MAX_VECTORS_TO_TEST = 15;
    static constexpr int UNREDUCED = -1;

    class PotentialPivotRef;

    class PotentialPivot {
    protected:
        int j; // best pivot col, or basic var if inactive
    public:
        double Sprime; // entropy of best pivot col / functionEntropy^|Mij|
        double R;
        double Rprime;
        typename std::multimap<double,int>::iterator bestPivotsEntry; // iterator into this row's bestPivots entry or value initialisaed if row is inactive

        PotentialPivot() : R(0.0) {}

        bool isScheduledForUpdate(std::multimap<double,int> &bestPivots) const { return bestPivotsEntry == bestPivots.end(); }
        void setScheduledForUpdate(std::multimap<double,int> &bestPivots) {
            bestPivots.erase(bestPivotsEntry);
            bestPivotsEntry = bestPivots.end();
        }
        bool isActive() const { return R != -1.0; }
        int basicVar() const {
            assert(!isActive());
            return j;
        }
        int bestPivotCol() const {
            assert(isActive());
            return j;
        }
        void setInactive(std::multimap<double,int> &bestPivots) {
            bestPivots.erase(bestPivotsEntry);
            bestPivotsEntry = typename std::multimap<double,int>::iterator();
            R = -1.0;
        }
        bool bestPivotIsKnown() const { return j >= 0; }
        void setBestPivotUnknown() { j = -1; }
        void setBestPivot(int j) {
            assert(isActive());
            this->j = j;
        }
//        const double &entropyChangeForecast() { return bestPivotsEntry->first; }
//        void recalculateEntropyChangeForecast(std::set<PotentialPivotRef> &bestPivots) {
//            assert(bestPivotsEntry == bestPivots.end());
//            bestPivots.emplace(calculateEntropyChangeForecast(), );
//            entropyChange = ;
//            std::cout << "Updated entropyChange to " << entropyChange << std::endl;
//            auto ins = bestPivots.insert(std::move(node));
//            assert(bestPivotsEntry == ins.position);
//            auto it = bestPivotsEntry;
//            it--;
//            assert(bestPivotsEntry == bestPivots.begin() || it->pivot.entropyChange < entropyChange);
//            for(const auto &entry: bestPivots) {
//                std::cout << entry.pivot.entropyChangeForecast() << " ";
//            }
//            std::cout << std::endl;
//        }
//        void initEntropyChangeForecast(std::set<PotentialPivotRef> &bestPivots) {
//            entropyChange = calculateEntropyChangeForecast();
//            bestPivotsEntry = bestPivots.emplace(*this).first;
//        }

        double calculateEntropyChangeForecast() const { return Sprime*(Rprime - Sprime) - R; }

        int basicColIndex() const { return isActive()? -1: j; }


    };

    class Row: public SortedSparseVec<COEFF> {
    public:
        double funcEntropy; // entropy associated with this row's function

        SortedSparseVec<COEFF> &coefficients() { return *this; }
        double functionEntropy(int Mij) const {
            return pow(funcEntropy, abs(Mij)); // TODO: Justify this!
        }
    };

//    class PotentialPivotRef {
//    public:
//        PotentialPivot &pivot;
//
//        PotentialPivotRef(PotentialPivot &pivot): pivot(pivot) { }
//
//        bool operator <(const PotentialPivotRef &other) const {
//            if(pivot.entropyChangeForecast() < other.pivot.entropyChangeForecast()) return true;
//            if(pivot.entropyChangeForecast() > other.pivot.entropyChangeForecast()) return false;
//            return &pivot < &other.pivot;
//        }
//    };

    class Column: public SortedSparseVec<COEFF> {
    public:
        double colEntropy;
        double funcEntropy;

        void setBasic() { colEntropy = -1.0; }
        bool isBasic() const { return signbit(colEntropy); }
    };

    std::vector<Column>     cols;   // tableau stored as col-major
    std::vector<Row>        rows;   // tableau stored as row-major (double storage is fastest way)
    std::multimap<double,int>   bestPivots; // entropy change forecst to row index
    std::vector<PotentialPivot> potentialPivots; // best pivot on each constraint row
    std::vector<COEFF>      F;     // constant in linear equations by row (see intro above)


//    template<class DOMAIN>
//    TableauEntropyMaximiser(const ConstrainedFactorisedDistribution<DOMAIN,COEFF> &distribution, double constraintKappa):
//            TableauEntropyMaximiser(calculateEntropiesByVarIndex(distribution), distribution.constraints, constraintKappa)
//            { }

//    template<class SAMPLER>
//    TableauEntropyMaximiser(const SAMPLER &priorSampler, const EqualityConstraints<COEFF> &constraints, double constraintKappa):
//            TableauEntropyMaximiser(calculateEntropiesByVarIndex(priorSampler), constraints, constraintKappa)
//    { }


    TableauEntropyMaximiser(const std::vector<double> &entropiesByVarIndex, const EqualityConstraints<COEFF> &constraints, double constraintKappa);

    void factorise();
    std::vector<SparseVec<COEFF>> getBasisVectors();
    EqualityConstraints<COEFF> getFactorisedConstraints();
//    template<class DOMAIN>
//    Basis<DOMAIN> getBasis();
    SparseVec<COEFF> getOrigin();

//    template<typename VEC>
//    void snapToSubspace(VEC &X) const;


//    std::pair<int,int> findBestPivot();

    void pivot(int i,int j);

    void sanityCheck();
protected:

    std::vector<int> prePivotUpdates(int i);
    void postPivotUpdates();

//    template<class DOMAIN>
//    static std::vector<double> calculateEntropiesByVarIndex(const ConstrainedFactorisedDistribution<DOMAIN,COEFF> &dist);

//    template<class SAMPLER>
//    static std::vector<double> calculateEntropiesByVarIndex(const SAMPLER &priorSampler);

    //    int sparsestColInRow(int i);
//    int sparsestRowInCol(int j);
//
//    void addRowSparsityEntry(int i);
//    void removeRowSparsityEntry(int i);
//    void addColSparsityEntry(int j);
//    void removeColSparsityEntry(int j);
//
//    void setColBasic(int j, int i);
//    void inactivateRow(int i);
//
//    double meanColumnL0Norm() const;
//    double meanColumnL1Norm() const;

    void postPivotUpdates(int i, int j, const std::vector<int> &affectedRows);

    double calculateR(int i);

    double calculateRPrime(int i);

    std::pair<int, double> findBestPivotOnRow(int i);

    double calculateColEntropy(int j);
};


// entropiesByVarIndex should be E_j = \sum_X (P(X)P(X+1_j))^0.5 = \sum_X P(X)(P(X+1_j)/P(X))^0.5
// where 1_j is a vector with a one in index j and zero elsewhere.
// constraints are the linear constraints
template<class COEFF>
TableauEntropyMaximiser<COEFF>::TableauEntropyMaximiser(const std::vector<double> &entropiesByVarIndex, const EqualityConstraints<COEFF> &constraints, double constraintKappa):
    cols(entropiesByVarIndex.size()),
    rows(constraints.size()),
    potentialPivots(constraints.size()),
    F(constraints.size())
{
    // initialise row entries
    for(int i=0; i < constraints.size(); ++i) {
        rows[i].coefficients() = constraints[i].coefficients;
        rows[i].funcEntropy = exp(-0.5*constraintKappa); // Entropy of F(x) = Ae^-k|x|
        F[i] = constraints[i].constant;
    }

    // initialise col entries
    for(int i=0; i<rows.size(); ++i) {
        for(const auto &[j,Mij]: rows[i]) {
            cols[j].emplace_back(i,Mij);
        }
    }
    for(int j=0; j<cols.size(); ++j) cols[j].sort();

    // initialise col entropies
    std::cout << "Col entropies" << std::endl;
    for(int j=0; j<cols.size(); ++j) {
        cols[j].funcEntropy = entropiesByVarIndex[j];
        cols[j].colEntropy = calculateColEntropy(j);
        std::cout << cols[j].colEntropy << " ";
    }
    std::cout << std::endl;

    // initialise potential pivots
    for(int i=0; i < constraints.size(); ++i) {
        potentialPivots[i].R = calculateR(i);
        potentialPivots[i].Rprime = calculateRPrime(i);
        std::pair<int,double> piv = findBestPivotOnRow(i);
        potentialPivots[i].setBestPivot(piv.first);
        potentialPivots[i].Sprime = piv.second;
        potentialPivots[i].bestPivotsEntry = bestPivots.emplace(potentialPivots[i].calculateEntropyChangeForecast(), i);
    }
    debug(sanityCheck());
}


// Use the zero vector to fit exponential approximation
//template<class COEFF>
//template<class DOMAIN>
//std::vector<double> TableauEntropyMaximiser<COEFF>::calculateEntropiesByVarIndex(const ConstrainedFactorisedDistribution<DOMAIN,COEFF> &dist) {
//    std::vector<double> entropiesByVarIndex
//}

// Use the zero vector to fit exponential approximation
//template<class COEFF>
//template<class SAMPLER>
//std::vector<double> TableauEntropyMaximiser<COEFF>::calculateEntropiesByVarIndex(const SAMPLER &priorSampler) {
//    std::cout << "Calculating dimension entropies" << std::endl;
//
//    int domainSize = priorSampler().size();
//    std::vector<int> sampleCounts(domainSize,0);
//    int targetMeanSampleCount = 100;
//    int targetTotalCount = domainSize*targetMeanSampleCount;
//    int nSamples = 0;
//    while(targetTotalCount > 0) {
//        const auto &sample = priorSampler();
//        for(int i=0; i<domainSize; ++i) {
//            targetTotalCount -= sample[i];
//            sampleCounts[i] += sample[i];
//        }
//        ++nSamples;
//    }
//
//    // fit exponential distribution to samples
//    std::vector<double> entropiesByVarIndex(domainSize);
//    for(int i=0; i<domainSize; ++i) {
//        double expki = sampleCounts[i]*1.0/(sampleCounts[i] + nSamples); // e^k of fitted exponential distribution
//        entropiesByVarIndex[i] = sqrt(expki); // entropy at Delta = 1
//    }
//
//    std::cout << entropiesByVarIndex << std::endl;
//    return entropiesByVarIndex;
//}


template<class COEFF>
void TableauEntropyMaximiser<COEFF>::factorise() {
    while(!bestPivots.empty()) {
        int bestPivotRow = bestPivots.rbegin()->second;
        int bestPivotCol = potentialPivots[bestPivotRow].bestPivotCol();
//        std::cout << "pivoting on " << bestPivotRow << ", " << bestPivotCol << " " << "colSize = " << cols[bestPivotCol].size() << " entropy change forecast = " << bestPivots.rbegin()->first << std::endl;
        pivot(bestPivotRow, bestPivotCol);
//        std::cout << *this << std::endl;
    }
}



template<class COEFF>
void TableauEntropyMaximiser<COEFF>::pivot(int i, int j) {
    Row &pivotRow = rows[i];
    Column &pivotCol = cols[j];
    COEFF minusMij = -pivotCol.get(i);
    assert(abs(minusMij) == 1);

    std::vector<int> affectedRows = prePivotUpdates(i);

    // pivot rows
    for(const std::pair<int,COEFF> &pivotColEntry: pivotCol) {
        const int &k = pivotColEntry.first;
        if(k != i) {
            COEFF multiplier = pivotColEntry.second/minusMij; // leave pivot element unchanged
            rows[k].weightedPlusAssign(multiplier, pivotRow);
            F[k] += multiplier * F[i];
        }
    }

    // pivot columns
    pivotCol.set(i,0); // leave pivot row unchanged (sort this out later)
    for(const std::pair<int,COEFF> &pivotRowEntry: pivotRow) {
        const int &l = pivotRowEntry.first;
        if(l != j) {
            cols[l].weightedPlusAssign(pivotRowEntry.second / minusMij, pivotCol);
        }
    }
    pivotCol.clear(); // set pivot column to pivot element only, with original value
    pivotCol.set(i, -minusMij);
    pivotCol.setBasic();

    postPivotUpdates(i, j, affectedRows);

    debug(sanityCheck());
}


// performs updates or row states before a pivot by
// subtracting columns whose contrib will change from R_i and R'_i
// removing rows whose entropy will change from bestPivots
// and recording which rows must be updated after the pivot
// Returns a vector of rows whose scores are affected by a pivot on row i
template<class COEFF>
std::vector<int> TableauEntropyMaximiser<COEFF>::prePivotUpdates(int i) {
    std::vector<int> affectedRows;
    potentialPivots[i].setInactive(bestPivots);
    affectedRows.reserve(rows[i].sparseSize() * cols[rows[i].front().first].sparseSize()); // approximate size
    for(const auto &[j,Mij]: rows[i]) {
        for(const auto &[k,Mkj]: cols[j]) {
            PotentialPivot &potentialPivotk = potentialPivots[k];
            if (potentialPivotk.isActive()) {
                potentialPivotk.R -= cols[j].colEntropy;
                potentialPivotk.Rprime -= cols[j].colEntropy / rows[k].functionEntropy(Mkj);
                if (!potentialPivotk.isScheduledForUpdate(bestPivots)) {
                    potentialPivotk.setScheduledForUpdate(bestPivots);
                    affectedRows.push_back(k);
                }
                if (potentialPivotk.bestPivotCol() == j) { // this col is best pivot for row k
                    potentialPivotk.setBestPivotUnknown();
                }
            }
        }
    }
    return affectedRows;
}


// recalculates all changed variables after a pivot. i.e.
//   the column entropies
//   the R and R' of each affected row
//   the best pivot point for each row
//   the row entropy change forecasts
//   the bestPivot entries
//
template<class COEFF>
void TableauEntropyMaximiser<COEFF>::postPivotUpdates(int i, int j, const std::vector<int> &affectedRows) {
    double R = 0.0;
    double Rprime = 0.0;
    double bestSij = -INFINITY;

    rows[i].funcEntropy = cols[j].funcEntropy; // pivot row now represents basic column

    // recalculate affected colEntropies
    for(const auto &[l,Mij]: rows[i]) {
        if(l != j) cols[l].colEntropy = calculateColEntropy(l);
    }

    // update best row pivots, R and R'
    for(const auto &[j,Mij]: rows[i]) {
        for(const auto &[k,Mkj]: cols[j]) {
            PotentialPivot &potentialPivotk = potentialPivots[k];
            if(potentialPivotk.isActive()) {
                potentialPivotk.R += cols[j].colEntropy;
                double Sprime = cols[j].colEntropy / rows[k].functionEntropy(Mkj);
                potentialPivotk.Rprime += Sprime;
                if (potentialPivotk.bestPivotIsKnown() && Sprime > potentialPivotk.Sprime) {
                    potentialPivotk.Sprime = Sprime;
                    potentialPivotk.setBestPivot(j);
                }
            }
        }
    }

    // now recalculate row entropies
    for(int k: affectedRows) {
        assert(k < potentialPivots.size());
        PotentialPivot &affectedRow = potentialPivots[k];
        if(!affectedRow.bestPivotIsKnown()) {
            affectedRow.Sprime = -INFINITY;
            for(const auto &[j,Mij]: rows[k]) {
                double Sprimej = cols[j].colEntropy / rows[k].functionEntropy(Mij);
                if(Sprimej > affectedRow.Sprime) {
                    affectedRow.Sprime = Sprimej;
                    affectedRow.setBestPivot(j);
                }
            }
        }
        assert(affectedRow.bestPivotIsKnown());
        assert(affectedRow.bestPivotCol() >=0);
//        std::cout << "updating entropy change for row " << k << std::endl;
        assert(affectedRow.bestPivotsEntry == bestPivots.end());
        affectedRow.bestPivotsEntry = bestPivots.emplace(affectedRow.calculateEntropyChangeForecast(), k);
    }


//    int bestCol = -1;
//    for(const std::pair<int,COEFF> &rowEntry: rows[i]) {
//        double Sj = colEntropy[rowEntry.first];
//        double Sij = Sj/rowEntry.second;
//        R += Sj;
//        Rprime += Sij;
//        if(Sij > bestSij) {
//            bestSij = Sij;
//            bestCol = rowEntry.first;
//        }
//    }
////    auto bestPivotNode = bestPivots.extract(bestPivotByRow[i]);
//
//    RowState &thisRowState = rowState[i];
//    if(thisRowState.setIterator != bestPivots.end()) {
//        bestPivots.erase(thisRowState.setIterator);
//        thisRowState.setIterator = bestPivots.end();
//    }
//    thisRowState.entropyChangeForecast = (Rprime - bestSij)*bestSij - R;
//    thisRowState.bestPivotj = bestCol;
//    bestPivots.insert(RowBestPivotRef(thisRowState));
}


template<class COEFF>
void TableauEntropyMaximiser<COEFF>::sanityCheck() {
    // check row/col are consistent
    for(int i=0; i<rows.size(); ++i) {
        assert(std::adjacent_find(rows[i].begin(), rows[i].end()) == rows[i].end());
        for(const auto &[j,Mij]: rows[i]) {
            assert(cols[j][i] == Mij);
        }
    }
    int nColEntries = 0;
    for(const auto &col: cols) { nColEntries += col.sparseSize(); }
    int nRowEntries = 0;
    for(const auto &row: rows) { nRowEntries += row.sparseSize(); }
    assert(nColEntries == nRowEntries);

    // check Potential pivot data is consistent
    int nActive = 0;
    for(int i=0; i < potentialPivots.size(); ++i) {
        if(potentialPivots[i].isActive()) {
            assert(fabs(potentialPivots[i].R - calculateR(i)) < 1e-8);
            assert(fabs(potentialPivots[i].Rprime - calculateRPrime(i)) < 1e-8);
            auto [bestPivj, bestSprime] = findBestPivotOnRow(i);
            assert(potentialPivots[i].bestPivotCol() == bestPivj);
            assert(potentialPivots[i].Sprime == bestSprime);
            assert(potentialPivots[i].bestPivotsEntry->second == i);
            assert(potentialPivots[i].bestPivotsEntry->first == potentialPivots[i].calculateEntropyChangeForecast());
            ++nActive;
        } else {
            assert(cols[potentialPivots[i].basicColIndex()].sparseSize() == 1);
        }
    }
    assert(nActive == bestPivots.size());

    // check colEntropies are consistent
    for(int j=0; j<cols.size(); ++j) {
        assert(cols[j].funcEntropy <= 1.0 && cols[j].funcEntropy >= 0.0);
        if(!cols[j].isBasic()) {
            assert(cols[j].colEntropy <= 1.0 && cols[j].colEntropy >= 0.0);
            assert(fabs(cols[j].colEntropy - calculateColEntropy(j)) < 1e-8);
        }
    }
//    std::cout << "passed sanity check" << std::endl;
}

// R_i = \sum_{j, Mij!=0} S_j
template<class COEFF>
double TableauEntropyMaximiser<COEFF>::calculateR(int i) {
    double R = 0.0;
    for(const auto &[j,Mij]: rows[i]) R += cols[j].colEntropy;
    return R;
}

// R'_i = \sum_{j, Mij!=0} S_j/F_i^|M_ij|
template<class COEFF>
double TableauEntropyMaximiser<COEFF>::calculateRPrime(int i) {
    double Rprime = 0.0;
    for(const auto &[j,Mij]: rows[i]) Rprime += cols[j].colEntropy / rows[i].functionEntropy(Mij);
    return Rprime;
}

template<class COEFF>
std::pair<int,double> TableauEntropyMaximiser<COEFF>::findBestPivotOnRow(int i) {
    int bestPiv = -1;
    double bestSprime = -INFINITY;
    for(const auto &[j,Mij]: rows[i]) {
        double Sprime = cols[j].colEntropy / rows[i].functionEntropy(Mij);
        if(Sprime > bestSprime) {
            bestSprime = Sprime;
            bestPiv = j;
        }
    }
    return {bestPiv, bestSprime};
}

template<class COEFF>
double TableauEntropyMaximiser<COEFF>::calculateColEntropy(int j) {
    assert(!cols[j].isBasic());
    double colEntropy = cols[j].funcEntropy;
    for (const auto &[k, Mkj]: cols[j]) {
        colEntropy *= rows[k].functionEntropy(Mkj);
    }
    return colEntropy;
}

template<class T>
std::vector<SparseVec<T>> TableauEntropyMaximiser<T>::getBasisVectors() {

    std::vector<SparseVec<T>> basisVectors;

    basisVectors.reserve(cols.size());
    for(int j=0; j<cols.size(); ++j) {
        if(!cols[j].isBasic()) {
            SparseVec<T> &newBasis = basisVectors.emplace_back();
            newBasis.insert(j, 1); // element from the identity (insert first to ensure canonical)
            for (auto & [i, Mij]: cols[j]) {
                int basicj = potentialPivots[i].basicVar();
                newBasis.insert(basicj, Mij / (-rows[i][basicj]));
            }
        }
    }
    return basisVectors;
}


// returns the value of the domain when all non-basic variables are zero
template<class T>
SparseVec<T> TableauEntropyMaximiser<T>::getOrigin() {
    SparseVec<T> origin;
    for(int row =0; row < potentialPivots.size(); ++row) {
        int basicVar = potentialPivots[row].basicVar();
        origin.insert(basicVar, F[row]/cols[basicVar][row]);
    }
    return origin;
}


template<class T>
std::ostream &operator <<(std::ostream &out, const TableauEntropyMaximiser<T> &tableau) {
    // print tableau
    for(int i=0; i<tableau.F.size(); ++i) {
        out << "x(" << tableau.potentialPivots[i].basicColIndex() << ")";
        out << "\t=\t";
        for(int j=0; j<tableau.cols.size(); ++j) {
            out << tableau.rows[i][j] << "\t";
        }
        out << "+\t" << tableau.F[i] << std::endl;
    }

    for(int i=tableau.F.size(); i<tableau.rows.size(); ++i) {
        out << "F_" << i-tableau.F.size();
        out << "  \t\t";
        for(int j=0; j<tableau.cols.size(); ++j) {
            out << tableau.rows[i][j] << "\t";
        }
        out << std::endl;
    }
    // mark basic columns
    out << "\t   \t \t";
    for(int j=0; j<tableau.cols.size(); ++j) {
        out << (tableau.cols[j].isBasic()?"B\t":"-\t");
    }
    out << std::endl;
    return out;
}

// returns (column index, change in entropy)
//template<class COEFF>
//std::pair<int,double> TableauEntropyMaximiser<COEFF>::getRowScore(int i) {
//    double R = 0.0;
//    double Rprime = 0.0;
//    double bestSij = -INFINITY;
//    int bestCol = -1;
//    for(const std::pair<int,COEFF> &rowEntry: rows[i]) {
//        double Sj = colEntropy[rowEntry.first];
//        double Sij = Sj/rowEntry.second;
//        R += Sj;
//        Rprime += Sij;
//        if(Sij > bestSij) {
//            bestSij = Sij;
//            bestCol = rowEntry.first;
//        }
//    }
//    auto bestPivotNode = bestPivots.extract(row[i]);
//    bestPivotNode.value().entropyChangeForecast = (Rprime - bestSij)*bestSij - R;
//    bestPivotNode.value().bestPivotj = bestCol;
//    bestPivots.insert(std::move(bestPivotNode));
//}

#endif //ABMCMC_TABLEAUENTROPYMAXIMISER_H
