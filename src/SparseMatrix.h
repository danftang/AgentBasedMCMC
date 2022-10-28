//
// Created by daniel on 17/10/22.
//

#ifndef ABMCMC_SPARSEMATRIX_H
#define ABMCMC_SPARSEMATRIX_H

#include <map>
#include <list>
#include <vector>
#include <set>
#include "include/SparseVec.h"

template<class T>
class SparseMatrix {
public:

//    class ColumnEntry;
//
//    class RowEntry {
//        T value;
//        std::list<ColumnEntry> &column;
//        std::map<int,RowEntry> &row;
//        typename std::list<ColumnEntry>::iterator ColumnEntryIt;
//        typename std::map<int,RowEntry>::iterator RowEntryIt;
//    };

    class Entry {
    public:
        const int i; // row index
        const int j; // col index
        T value;
    };

    class EntryRef: public std::list<Entry>::iterator {
    public:
        bool operator <(const EntryRef &other) const { return this->j < other->j; }
    };

    class Column: public std::set<EntryRef> { // set of row indices of non-zero entries, +ve is equality constraint, -ve is factor dependency
    public:
    };


    class Row: public std::list<Entry> { // map from non-zero col index to element value
    public:

        Row(const SparseVec<T> &row);
        Row() {}

        Row &operator *=(const T &c);

        // dot product
        template<typename OTHER, typename E = decltype(std::declval<T &>() = std::declval<T>() * std::declval<OTHER>()[0])>
        T operator *(const OTHER &other) const {
            T dotProd = 0;
            for(const auto &element: *this) {
                dotProd += element.second * other[element.first];
            }
            return dotProd;
        }
    };

    std::vector<Row>    rows;
    std::vector<Column> cols;

    Row & operator [](int i) {
        assert(i < rows.size());
        return rows[i];
    }

    T & at(int i, int j) {
        assert(i < rows.size() && j < cols.size());
        return rows[i].at(j).value;
    }


    void pivot(int pi, int pj);


};

template<class T>
void SparseMatrix<T>::pivot(int pi, int pj) {
    const T &Mpipj = rows[pi].at(pj);

    rows[pi] *= -1.0/Mpipj;     // divide row through to make the pivot point -1

    for (int i: cols[pj]) {
        if (i != pi) {
            if (i >= 0) { // row is equality constraint, so do linear pivot
                T Mipj = rows[i].at(pj);
                for (const auto &[j, Mpij]: rows[pi]) {
                    if (j != pj) {
                        auto MijIt = rows[i].find(j);
                        if (MijIt == rows[i].end()) { // fill-in
                            rows[i][j] = Mipj * Mpij;
                            cols[j].insert(i);
                        } else {
                            T newMij = (MijIt->second += Mipj * Mpij);
                            if (newMij == 0) { // drop-out
                                cols[j].erase(i);
                                rows[i].erase(j);
                            }
                        }
                    } else {
                        // drop out without testing as we're on the pivot col
                        // don't erase col entry yet as we're iterating through it
                        rows[i].erase(j);
                    }
                }
            } else { // row is factor dependency so do union pivot
                for (const auto &[j, Mpij]: rows[pi]) cols[j].insert(i);
            }
        }
    }

//    // re-add previously removed sparsity entries
//    for(int i : cols[pj]) {
//        if(i != pi) {
//            rows[i].erase(pj);
//            if(rows[i].isActive) addRowSparsityEntry(i);
//        }
//    }
    setColBasic(pj, pi); // removes entries from col and records in basis
    // add new column sparsity entries
    for(const auto [j, Mpij]: rows[pi]) {
        if(j != pj) {
//            cols[j].erase(pi);  // cols[] only stores non-reduced rows (triangular factorisation)
            addColSparsityEntry(j);
        }
    }
    inactivateRow(pi); // sets inactive and removes row sparsity entry
}


#endif //ABMCMC_SPARSEMATRIX_H
