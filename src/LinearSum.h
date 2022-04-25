// Represents a sum of terms sum_i c_ix_i
//
// Created by daniel on 26/04/2021.
//

#ifndef ABMCMC_LINEARSUM_H
#define ABMCMC_LINEARSUM_H

#include <map>
#include <iomanip>
#include "X.h"
#include "SparseVec.h"

template<class T>
class LinearSum {
public:
    std::map<int,T> coefficients;

    LinearSum() {}

    LinearSum(std::initializer_list<std::pair<const int,T>> initializerList) : coefficients(initializerList) { }

    LinearSum(LinearSum<T> &&rvalue): coefficients(std::move(rvalue.coefficients)) { }

    LinearSum(const LinearSum<T> &lvalue) : coefficients(lvalue.coefficients) { }

    T &operator [](int index) {
        typename std::map<int,T>::iterator it = coefficients.find(index);
        if(it == coefficients.end()) return coefficients[index] = 0; // create new entry and init to zero if non existant
        return it->second;
    }

    LinearSum<T> &operator +=(const LinearSum<T> &other) {
        for(auto [varId, coeff]: other.coefficients) {
            (*this)[varId] += coeff;
        }
        return *this;
    }

    LinearSum<T> &operator -=(const LinearSum<T> &other) {
        for(auto [varId, coeff]: other.coefficients) {
            (*this)[varId] -= coeff;
        }
        return *this;
    }

    LinearSum<T> &operator *=(T multiplier) {
        for(auto it= coefficients.begin(); it != coefficients.end(); ++it) {
            it->second *= multiplier;
        }
        return *this;
    }

    LinearSum<T> operator +(const LinearSum<T> &rhs) {
        LinearSum sum(*this);
        sum += rhs;
        return sum;
    }

    LinearSum<T> operator +(LinearSum<T> &&rhs) {
        LinearSum sum(rhs);
        sum += *this;
        return sum;
    }

    LinearSum<T> operator -(const LinearSum<T> &rhs) {
        LinearSum lhs(*this);
        lhs -= rhs;
        return lhs;
    }

    LinearSum<T> operator -(LinearSum<T> &&rhs) {
        rhs *= -1;
        return rhs + *this;
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
inline LinearSum<T> operator *(T coefficient, const X &variable) {
    return LinearSum<T>({ {variable.id, coefficient} });
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
