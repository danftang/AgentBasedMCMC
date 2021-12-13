//
// Created by daniel on 10/12/2021.
//

#include "TableauNormMinimiser.h"


TableauNormMinimiser::TableauNormMinimiser(glp::Problem &problem) {
    rows.reserve(problem.nConstraints());
    for(int i=1; i<=problem.nConstraints(); ++i) {
        bool isActive = (problem.getRowType(i) == GLP_FX);
        rows.push_back(Row(problem.getMatRow(i), isActive));
        rowsBySparsity[rows.back().sparseSize].push_front(i-1);
        rows.back().sparsityEntry = rowsBySparsity[rows.back().sparseSize].begin();
    }
    cols.reserve(problem.nVars());
    for(int j=1; j<=problem.nVars(); ++j) {
        cols.push_back(Column(problem.getMatCol(j)));
        colsBySparsity[cols.back().sparseSize].push_front(j-1);
        cols.back().sparsityEntry = colsBySparsity[cols.back().sparseSize].begin();
    }
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
    cols[pj].setBasic(pi, colsBySparsity);
    for(const auto [j, Mpij]: rows[pi]) updateColSparsity(j);
}


void TableauNormMinimiser::updateRowSparsity(int i) {
    Row &row = rows[i];
    if(row.sparseSize != row.size()) {
        rowsBySparsity[row.sparseSize].erase(row.sparsityEntry);
        row.sparseSize = row.size();
        rowsBySparsity[row.sparseSize].push_front(i);
        row.sparsityEntry = rowsBySparsity[row.sparseSize].begin();
    }
}


void TableauNormMinimiser::updateColSparsity(int j) {
    Column &col = cols[j];
    if(col.sparseSize != col.size()) {
        colsBySparsity[col.sparseSize].erase(col.sparsityEntry);
        col.sparseSize = col.size();
        colsBySparsity[col.sparseSize].push_front(j);
        col.sparsityEntry = rowsBySparsity[col.sparseSize].begin();
    }
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

void TableauNormMinimiser::Column::setBasic(int pi, std::vector<std::list<int>> &colsBySparsity) {
    isBasic = true;
    colsBySparsity[sparseSize].erase(sparsityEntry);
    sparsityEntry = colsBySparsity[sparseSize].end();
    clear();
    insert(pi);
    sparseSize = 1;
}
