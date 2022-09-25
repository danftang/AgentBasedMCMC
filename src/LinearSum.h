// Represents a sum of terms sum_i c_ix_i
//
// Created by daniel on 26/04/2021.
//

#ifndef ABMCMC_LINEARSUM_H
#define ABMCMC_LINEARSUM_H

#include <map>
#include <iomanip>
#include "X.h"
#include "include/SparseVec.h"

template<class T>
class LinearSum {
public:
    std::map<int,T> coefficients;

    LinearSum()=default;

    LinearSum(std::initializer_list<std::pair<const int,T>> initializerList) : coefficients(initializerList) { }

    LinearSum(LinearSum<T> &&rvalue): coefficients(std::move(rvalue.coefficients)) { }

    LinearSum(const LinearSum<T> &lvalue) : coefficients(lvalue.coefficients) { }


    T &operator [](int index) {
        typename std::map<int,T>::iterator it = coefficients.find(index);
        if(it == coefficients.end()) return coefficients[index] = 0; // create new entry and init to zero if non existant
        return it->second;
    }

    template<class OTHER>
    LinearSum<T> &operator +=(const LinearSum<OTHER> &other) {
        for(auto [varId, coeff]: other.coefficients) {
            (*this)[varId] += coeff;
        }
        return *this;
    }

    template<class OTHER>
    LinearSum<T> &operator -=(const LinearSum<OTHER> &other) {
        for(auto [varId, coeff]: other.coefficients) {
            (*this)[varId] -= coeff;
        }
        return *this;
    }

    LinearSum<T> &operator +=(const X &rhs) {
        auto ins = coefficients.insert({rhs.id,1});
        if(!ins.second) ins.first->second += 1;
        return *this;
    }

    LinearSum<T> &operator -=(const X &rhs) {
        auto ins = coefficients.insert({rhs.id,-1});
        if(!ins.second) ins.first->second -= 1;
        return *this;
    }

    template<class MULT>
    LinearSum<T> &operator *=(MULT multiplier) {
        for(auto it= coefficients.begin(); it != coefficients.end(); ++it) {
            it->second *= multiplier;
        }
        return *this;
    }

    template<class OTHER, class RESULT = decltype(std::declval<T>() + std::declval<OTHER>())>
    LinearSum<RESULT> operator +(const LinearSum<OTHER> &rhs) {
        LinearSum<RESULT> sum;
        for(const auto &entry: coefficients) sum.coefficients[entry.first] = entry.second;
        for(const auto &entry: rhs.coefficients) {
            auto insert = sum.coefficients.insert(entry);
            if(!insert.second) insert.first->second += entry.second;
        }
        return sum;
    }

    template<class OTHER, class RESULT = decltype(std::declval<T>() + std::declval<OTHER>())>
    LinearSum<RESULT> operator -(const LinearSum<OTHER> &rhs) {
        LinearSum<RESULT> difference;
        for(const auto &entry: coefficients) difference.coefficients[entry.first] = entry.second;
        for(const auto &entry: rhs.coefficients) {
            auto insert = difference.coefficients.insert(std::pair(entry.first, 0));
            insert.first->second += entry.second;
        }
        return difference;
    }


    LinearSum<T> operator +(const X &rhs) && {
        (*this) += rhs;
        return std::move(*this);
    }


    LinearSum<T> operator -(const X &rhs) && {
        (*this) -= rhs;
        return std::move(*this);
    }

    LinearSum<T> operator +(const LinearSum<T> &rhs) && {
        (*this) += rhs;
        return std::move(*this);
    }

    LinearSum<T> operator -(const LinearSum<T> &rhs) && {
        (*this) -= rhs;
        return std::move(*this);
    }

    SparseVec<T> toSparseVec() const {
        SparseVec<T> v(coefficients.size());
        int i=0;
        for(auto [key,val] : coefficients) {
            v.indices[i] = key;
            v.values[i++] = val;
        }
        return v;
    }

};

template<class T>
inline LinearSum<T> operator *(const T &coefficient, const X &variable) {
    return LinearSum<T>({ {variable.id, coefficient} });
}

// If coefficient type is uncertain, go for the lowest precedence type
// w.r.t. arithmetic operations
inline LinearSum<char> operator +(const X &lhsX, const X &rhsX) {
    return LinearSum<char>({ {lhsX.id, 1}, {rhsX.id, 1} });
}

inline LinearSum<char> operator -(const X &lhsX, const X &rhsX) {
    return LinearSum<char>({ {lhsX.id, 1}, {rhsX.id, -1} });
}



template<class T>
std::ostream &operator <<(std::ostream &out, const LinearSum<T> &sum) {
    bool first = true;
    for (auto [varId,coeff] : sum.coefficients) {
        if(first) first = false; else out << " +\t";
        out << std::setw(8) << coeff << "X" << varId;
    }
    return out;
}

#endif //GLPKPP_LAZYLINEAREXPRESSION_H
