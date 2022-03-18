//
// Created by daniel on 10/12/2021.
//

#include <iostream>
#include <assert.h>
#include <cmath>
#include "StlStream.h"
#include "TableauNormMinimiser.h"

// Turns the supplied convex polyhedron into a set of equality constraints by introducing auxiliary
// variables where the upper and lower bounds are not equal. So if we have
// L <= CX <= U
// then
// (-L M)(1 X A)^T = 0, 0 <= A <= U-L
// where the zero'th elelemet of (1 X A) is 1.
TableauNormMinimiser::TableauNormMinimiser(const ConvexPolyhedron &problem)
{
    basis.reserve(problem.size());
    F.reserve(problem.size());
    rows.reserve(problem.size());
    Ha.reserve(problem.size());
    Ha.push_back(0.0);
    int maxCol = 0; // maximum column id seen so far
    nAuxiliaryVars = 0;
    for(const Constraint &constraint: problem) {
        if(constraint.coefficients.sparseSize() == 1) { // rectilinear limits
            assert(constraint.lowerBound == 0);
            int xi = constraint.coefficients.indices[0];
            if(Hc.size() < xi + 1) Hc.resize(xi+1);
            Hc[xi] = constraint.upperBound/constraint.coefficients.values[0];
            if(xi > maxCol) maxCol = xi;
        } else { // normal tableau entry
            bool isActive = (constraint.upperBound == constraint.lowerBound);
            rows.push_back(Row(constraint.coefficients, isActive));
            if (isActive) {
                addRowSparsityEntry(rows.size()-1);
                basis.push_back(0); // fixed var
            } else {
                basis.push_back(-(++nAuxiliaryVars));
                Ha.push_back(constraint.upperBound - constraint.lowerBound);
            }
            F.push_back(-constraint.lowerBound); // lower bounds are all set to 0
            int rowMaxCol = rows.back().rbegin()->first;
            if(rowMaxCol > maxCol) maxCol = rowMaxCol;
        }
    }

    Hc.resize(maxCol+1);
    cols.resize(maxCol +1);
    for(int j=1; j<=maxCol; ++j) assert(Hc[j] != 0); // no rectilinear limit given for non-basic j or fixed at 0
    for(int i=0; i<rows.size(); ++i) {
        for(auto [j, v]: rows[i]) {
            cols[j].insert(i);
        }
    }
    for(int j=0; j<=maxCol; ++j) {
        addColSparsityEntry(j);
    }
    std::cout << "Initial mean L0 norm = " << meanColumnL0Norm() << std::endl;
    std::cout << "Initial L1 norm = " << meanColumnL1Norm() << std::endl;
}


void TableauNormMinimiser::findMinimalBasis() {
    while(!rowsBySparsity.empty()) {
        auto [i,j] = findMarkowitzPivot();
        pivot(i,j);
//        std::cout << minimalBasis << std::endl;
//        std::cout << rowsBySparsity << std::endl;
    }
    std::cout << "Reduced mean L0 norm = " << meanColumnL0Norm() << std::endl;
    std::cout << "Reduced L1 norm = " << meanColumnL1Norm() << std::endl;

}


std::pair<int,int> TableauNormMinimiser::findMarkowitzPivot() {
    int k=1;
    int lowestScore = INT32_MAX;
    int vectorsTested = 0;
    std::pair<int,int> bestPivot(-1,-1);
    while(k<colsBySparsity.size() && k<rowsBySparsity.size()) {
//        std::cout << "k = " << k << std::endl;
        for(int j: colsBySparsity[k]) {
            int bestRow = sparsestRowInCol(j);
            if(bestRow != -1) {
                int score = (rows[bestRow].size() - 1) * (cols[j].size() - 1);
                if (score < lowestScore) {
                    bestPivot = std::pair<int, int>(bestRow, j);
                    if (score <= (k - 1) * (k - 1)) return bestPivot;
                    lowestScore = score;
                }
                if (++vectorsTested >= maxVectorsToTest) return bestPivot;
            }
        }
        for(int i: rowsBySparsity[k]) {
            int bestCol = sparsestColInRow(i);
            if(bestCol != -1) {
                int score = (rows[i].size() - 1) * (cols[bestCol].size() - 1);
                if (score < lowestScore) {
                    bestPivot = std::pair<int, int>(i, bestCol);
                    if (score <= (k - 1) * (k - 1)) return bestPivot;
                    lowestScore = score;
                }
                if (++vectorsTested >= maxVectorsToTest) return bestPivot;
            }
        }
        ++k;
    }
    return bestPivot;
}


void TableauNormMinimiser::pivot(int pi, int pj) {
    assert(rows[pi].isActive);
//    std::cout << "Pivoting at " << pi << " " << pj << std::endl;
    double Mpipj = rows[pi].at(pj);

    // remove sparsity entries of all affected rows and columns
    for( auto [j, v] : rows[pi]) removeColSparsityEntry(j);

    rows[pi] *= -1.0/Mpipj;     // divide row through to make the pivot point -1
    F[pi] *= -1.0/Mpipj;

    for(int i : cols[pj]) {
        if(i != pi) {
            if(rows[i].isActive) removeRowSparsityEntry(i);
            double Mipj = rows[i].at(pj);
//            double MipjOverMpipj = rows[i].at(pj)/Mpipj;
            for (const auto[j, Mpij]: rows[pi]) {
                if(j != pj) {
                    auto MijIt = rows[i].find(j);
                    if (MijIt == rows[i].end()) { // fill-in
                        rows[i][j] = Mipj * Mpij;
                        cols[j].insert(i);
                    } else {
                        double newMij = (rows[i][j] += Mipj * Mpij);
                        if (newMij == 0.0) { // drop-out
                            cols[j].erase(i);
                            rows[i].erase(j);
                        }
                    }
                }
            }
            F[i] += Mipj * F[pi];
        }
    }

    for(int i : cols[pj]) {
        if(i != pi) {
            rows[i].erase(pj);
            if(rows[i].isActive) addRowSparsityEntry(i);
        }
    }
    setColBasic(pj, pi);
    for(const auto [j, Mpij]: rows[pi]) {
        if(j != pj) addColSparsityEntry(j);
    }
    inactivateRow(pi);
}


//void TableauNormMinimiser::updateRowSparsity(int i) {
//    Row &row = rows[i];
//    if(row.sparseSize != row.size()) {
//        removeRowSparsityEntry(i);
//        addRowSparsityEntry(i);
//    }
//}



void TableauNormMinimiser::addRowSparsityEntry(int i) {
    Row &row = rows[i];
    if(rowsBySparsity.size() <= row.size()) rowsBySparsity.resize(row.size()+1);
    rowsBySparsity[row.size()].push_front(i);
    row.sparsityEntry = rowsBySparsity[row.size()].begin();
}

void TableauNormMinimiser::removeRowSparsityEntry(int i) {
    Row &row = rows[i];
    assert(row.isActive);
    rowsBySparsity[row.size()].erase(row.sparsityEntry);
    while(rowsBySparsity.size() > 0 && rowsBySparsity.back().empty()) rowsBySparsity.pop_back();
}



//void TableauNormMinimiser::updateColSparsity(int j) {
//    Column &col = cols[j];
//    if(col.sparseSize != col.size()) {
//        removeColSparsityEntry(j);
//        addColSparsityEntry(j);
//    }
//}

void TableauNormMinimiser::addColSparsityEntry(int j) {
    Column &col = cols[j];
//    col.sparseSize = col.size();
    if(colsBySparsity.size() <= col.size()) colsBySparsity.resize(col.size()+1);
    colsBySparsity[col.size()].push_front(j);
    col.sparsityEntry = colsBySparsity[col.size()].begin();
}

void TableauNormMinimiser::removeColSparsityEntry(int j) {
    Column &col = cols[j];
    assert(!col.isBasic);
    colsBySparsity[col.size()].erase(col.sparsityEntry);
    col.sparsityEntry = colsBySparsity[col.size()].end();
    while(colsBySparsity.size() > 0 && colsBySparsity.back().empty()) colsBySparsity.pop_back();
}


TableauNormMinimiser::Row::Row(const SparseVec<double> &row, bool isActive): isActive(isActive) {
    for (int k = 0; k < row.sparseSize(); ++k) {
        (*this)[row.indices[k]] = row.values[k];
    }
}

TableauNormMinimiser::Row &TableauNormMinimiser::Row::operator*=(double c) {
    for(auto it = begin(); it != end(); ++it) it->second *= c;
    return *this;
}


//TableauNormMinimiser::Column::Column(const SparseVec<double> &col): isBasic(false) {
//    for (int k = 0; k < col.sparseSize(); ++k) {
//        insert(col.indices[k]);
//    }
//}


void TableauNormMinimiser::setColBasic(int j, int i) {
    Column &col = cols[j];
    col.isBasic = true;
    col.clear();
    col.insert(i);
    basis[i] = j;
}

void TableauNormMinimiser::inactivateRow(int i) {
    Row &row = rows[i];
    removeRowSparsityEntry(i);
    row.isActive = false;
}

int TableauNormMinimiser::sparsestColInRow(int i) {
    int minColSparsity = INT32_MAX;
    int minj = -1;
    for(const auto [j, Mij] : rows[i]) {
        int colSize = cols[j].size();
        if(colSize < minColSparsity && !cols[j].isBasic && fabs(Mij) == 1.0) {
            minColSparsity = colSize;
            minj = j;
        }
    }
//    std::cout << "Considering pivot from row " << i << " " << minj << std::endl;
    return minj;
}

int TableauNormMinimiser::sparsestRowInCol(int j) {
    int minRowSparsity = INT32_MAX;
    int mini = -1;
    for(int i : cols[j]) {
        int rowSize = rows[i].size();
        if(rowSize < minRowSparsity && rows[i].isActive && fabs(rows[i].at(j)) == 1.0) {
            minRowSparsity = rowSize;
            mini = i;
        }
    }
//    std::cout << "Considering pivot " << mini << " from col " << j << std::endl;
    return mini;
}

double TableauNormMinimiser::meanColumnL0Norm() const {
    double normSum = 0.0;
    int nNonBasic = 0;
    for(int j=0; j<cols.size(); ++j) {
        if(!cols[j].isBasic) {
            normSum += cols[j].size();
            ++nNonBasic;
        }
    }
    return normSum / nNonBasic;
}

double TableauNormMinimiser::meanColumnL1Norm() const {
    double normSum = 0.0;
    int nNonBasic = 0;
    for(int j=0; j<cols.size(); ++j) {
        if(!cols[j].isBasic) ++nNonBasic;
    }
    for(int i=0; i<rows.size(); ++i) {
        for(const auto [j, Mij] : rows[i]) {
            if(!cols[j].isBasic) normSum += fabs(Mij);
        }
    }
    return normSum / nNonBasic;
}

std::ostream &operator <<(std::ostream &out, const TableauNormMinimiser &tableau) {
    // show upper limits
    out << "    \t \t";
    for(int j=0; j<tableau.cols.size(); ++j) {
        out << tableau.Hc[j] << "\t";
    }
    out << std::endl;
    // print tableau
    for(int i=0; i<tableau.rows.size(); ++i) {
        out << "x(" << tableau.basis[i] << ")\t=\t";
        for(int j=0; j<tableau.cols.size(); ++j) {
            out << tableau(i,j) << "\t";
        }
        out << "+\t" << tableau.F[i];
        if(tableau.basis[i] < 0) out << " <= " << tableau.Ha[-tableau.basis[i]];
        out << std::endl;
    }
    // mark basic columns
    out << "    \t \t";
    for(int j=0; j<tableau.cols.size(); ++j) {
        out << (tableau.cols[j].isBasic?"B\t":"-\t");
    }
    out << std::endl;
    return out;
}