//
// Created by daniel on 10/12/2021.
//

#include "TableauNormMinimiser.h"


TableauNormMinimiser::TableauNormMinimiser(glp::Problem &problem): minimalBasis(problem.nConstraints()) {
    rows.reserve(problem.nConstraints());
    for(int i=1; i<=problem.nConstraints(); ++i) {
        bool isActive = (problem.getRowType(i) == GLP_FX);
        if(!isActive) minimalBasis[i-1] = -i;
        rows.push_back(Row(problem.getMatRow(i), isActive));
        addRowSparsityEntry(i-1);
    }
    cols.reserve(problem.nVars());
    for(int j=1; j<=problem.nVars(); ++j) {
        cols.push_back(Column(problem.getMatCol(j)));
        addColSparsityEntry(j-1);
    }
}

void TableauNormMinimiser::findMinimalBasis() {
    while(!rowsBySparsity.empty()) {
        auto [i,j] = findMarkowitzPivot();
        pivot(i,j);
    }
}


std::pair<int,int> TableauNormMinimiser::findMarkowitzPivot() {
    int k=1;
    int lowestScore = INT32_MAX;
    int pivotsTested = 0;
    std::pair<int,int> bestPivot;
    while(k<colsBySparsity.size() && k<rowsBySparsity.size()) {
        for(int j: colsBySparsity[k]) {
            int bestRow = sparsestRowInCol(j);
            if(bestRow != -1) {
                int score = (rows[bestRow].sparseSize - 1) * (cols[j].sparseSize - 1);
                if (score < lowestScore) {
                    bestPivot = std::pair<int, int>(bestRow, j);
                    if (score <= (k - 1) * (k - 1)) return bestPivot;
                    lowestScore = score;
                }
                if (++pivotsTested >= maxPivotsToTest) return bestPivot;
            }
        }
        for(int i: rowsBySparsity[k]) {
            int bestCol = sparsestColInRow(i);
            if(bestCol != -1) {
                int score = (rows[i].sparseSize - 1) * (cols[bestCol].sparseSize - 1);
                if (score < lowestScore) {
                    bestPivot = std::pair<int, int>(i, bestCol);
                    if (score <= (k - 1) * (k - 1)) return bestPivot;
                    lowestScore = score;
                }
                if (++pivotsTested >= maxPivotsToTest) return bestPivot;
            }
        }
        ++k;
    }
    return std::pair<int,int>(-1,-1); // no valid pivots exist
}


void TableauNormMinimiser::pivot(int pi, int pj) {
    double Mpipj = rows[pi].at(pj);
    for(int i : cols[pj]) {
        if(i != pi) {
            double Mipj = rows[i].at(pj);
            for (const auto[j, Mpij]: rows[pi]) {
                auto MijIt = rows[i].find(j);
                if (MijIt == rows[i].end()) { // fill-in
                    rows[i][j] = -Mipj * Mpij / Mpipj;
                    cols[j].insert(i);
                } else {
                    double newMij = (rows[i][j] -= Mipj * Mpij / Mpipj);
                    if (newMij == 0.0) { // drop-out
                        cols[j].erase(i);
                        rows[i].erase(j);
                    }
                }
            }
        }
    }

    for(int i : cols[pj]) {
        rows[i].erase(pj);
        updateRowSparsity(i);
    }
    setColBasic(pj, pi);
    for(const auto [j, Mpij]: rows[pi]) updateColSparsity(j);
    inactivateRow(pi);
}


void TableauNormMinimiser::updateRowSparsity(int i) {
    Row &row = rows[i];
    if(row.sparseSize != row.size()) {
        removeRowSparsityEntry(i);
        addRowSparsityEntry(i);
    }
}



void TableauNormMinimiser::addRowSparsityEntry(int i) {
    Row &row = rows[i];
    row.sparseSize = row.size();
    if(rowsBySparsity.size() <= row.sparseSize) rowsBySparsity.resize(row.sparseSize+1);
    rowsBySparsity[row.sparseSize].push_front(i);
    row.sparsityEntry = rowsBySparsity[row.sparseSize].begin();
}

void TableauNormMinimiser::removeRowSparsityEntry(int i) {
    Row &row = rows[i];
    rowsBySparsity[row.sparseSize].erase(row.sparsityEntry);
    while(rowsBySparsity.back().empty()) rowsBySparsity.pop_back();
}



void TableauNormMinimiser::updateColSparsity(int j) {
    Column &col = cols[j];
    if(col.sparseSize != col.size()) {
        removeColSparsityEntry(j);
        addColSparsityEntry(j);
    }
}

void TableauNormMinimiser::addColSparsityEntry(int j) {
    Column &col = cols[j];
    col.sparseSize = col.size();
    if(colsBySparsity.size() <= col.sparseSize) colsBySparsity.resize(col.sparseSize+1);
    colsBySparsity[col.sparseSize].push_front(j);
    col.sparsityEntry = colsBySparsity[col.sparseSize].begin();
}

void TableauNormMinimiser::removeColSparsityEntry(int j) {
    Column &col = cols[j];
    colsBySparsity[col.sparseSize].erase(col.sparsityEntry);
    while(colsBySparsity.back().empty()) colsBySparsity.pop_back();
}


TableauNormMinimiser::Row::Row(const glp::SparseVec &row, bool isActive): isActive(isActive) {
    for (int k = 0; k < row.sparseSize(); ++k) {
        (*this)[row.indices[k]] = row.values[k];
    }
    sparseSize = row.sparseSize();
}


TableauNormMinimiser::Column::Column(const glp::SparseVec &col): isBasic(false) {
    insert(col.indices.begin(), col.indices.end());
    sparseSize = col.sparseSize();
}


void TableauNormMinimiser::setColBasic(int j, int i) {
    Column &col = cols[j];
    col.isBasic = true;
    removeColSparsityEntry(j);
    col.clear();
    col.insert(i);
    col.sparseSize = 1;
    minimalBasis[i] = j;
}

void TableauNormMinimiser::inactivateRow(int i) {
    Row &row = rows[i];
    row.isActive = false;
    removeRowSparsityEntry(i);
}

int TableauNormMinimiser::sparsestColInRow(int i) {
    int minColSparsity = INT32_MAX;
    int minj = -1;
    for(const auto [j, Mij] : rows[i]) {
        int colSize = cols[j].sparseSize;
        if(colSize < minColSparsity && !cols[j].isBasic && fabs(Mij) == 1.0) {
            minColSparsity = colSize;
            minj = j;
        }
    }
    return minj;
}

int TableauNormMinimiser::sparsestRowInCol(int j) {
    int minRowSparsity = INT32_MAX;
    int mini = -1;
    for(int i : cols[j]) {
        int rowSize = rows[i].sparseSize;
        if(rowSize < minRowSparsity && rows[i].isActive && fabs(rows[i].at(j)) == 1.0) {
            minRowSparsity = rowSize;
            mini = i;
        }
    }
    return mini;
}