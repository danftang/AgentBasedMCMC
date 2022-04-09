
//
// Created by daniel on 26/04/2021.
//

#ifndef ABMCMC_CONSTRAINT_H
#define ABMCMC_CONSTRAINT_H

#include "SparseVec.h"
#include "LinearSum.h"

template<class COEFF>
class Constraint {
public:

    SparseVec<COEFF> coefficients;
    COEFF upperBound;
    COEFF lowerBound;

    Constraint(): Constraint(std::numeric_limits<COEFF>::lowest(), std::numeric_limits<COEFF>::max()) {}
    Constraint(COEFF lowerBound, COEFF upperBound);
    Constraint(COEFF lowerBound, SparseVec<COEFF> coefficients, COEFF upperBound);

    Constraint & operator <=(COEFF upperBound);

    bool isValidSolution(const std::vector<COEFF> &X) const {
        COEFF mx = coefficients * X;
        return (lowerBound <= mx) && (mx <= upperBound);
    }
};


template<class COEFF>
Constraint<COEFF> &Constraint<COEFF>::operator<=(COEFF upperBound) {
    if(upperBound < this->upperBound) this->upperBound = upperBound;
    return *this;
}


template<class COEFF>
Constraint<COEFF>::Constraint(COEFF lowerBound, COEFF upperBound):
        upperBound(upperBound),
        lowerBound(lowerBound) { }


template<class COEFF>
Constraint<COEFF>::Constraint(COEFF lowerBound, SparseVec<COEFF> sum, COEFF upperBound):
        coefficients(std::move(sum)),
        upperBound(upperBound),
        lowerBound(lowerBound) { }


template<class COEFF>
Constraint<COEFF> operator ==(const LinearSum<COEFF> &linExp, COEFF c) {
    return Constraint(c, linExp.toSparseVec(), c);
}

template<class COEFF>
Constraint<COEFF> operator ==(double c, const LinearSum<COEFF> &linExp) {
    return linExp == c;
}

template<class COEFF>
Constraint<COEFF> operator <=(const LinearSum<COEFF> &linExp, COEFF c) {
    return Constraint(std::numeric_limits<COEFF>::lowest(), linExp.toSparseVec(), c);
}

template<class COEFF>
Constraint<COEFF> operator >=(const LinearSum<COEFF> &linExp, COEFF c) {;
    return Constraint(c, linExp.toSparseVec(), std::numeric_limits<COEFF>::max());
}

template<class COEFF>
Constraint<COEFF> operator <=(COEFF c, const LinearSum<COEFF> &linExp) {
    return linExp >= c;
}

template<class COEFF>
Constraint<COEFF> operator >=(COEFF c, const LinearSum<COEFF> &linExp) {
    return linExp <= c;
}

template<class COEFF>
std::enable_if_t<!std::is_same_v<COEFF,int>, Constraint<COEFF>>
operator <=(const LinearSum<COEFF> &linExp, int c) {
    return Constraint(std::numeric_limits<COEFF>::lowest(), linExp.toSparseVec(), COEFF(c));
}

template<class COEFF>
std::enable_if_t<!std::is_same_v<COEFF,int>, Constraint<COEFF>>
operator >=(const LinearSum<COEFF> &linExp, int c) {;
    return Constraint(COEFF(c), linExp.toSparseVec(), std::numeric_limits<COEFF>::max());
}

template<class COEFF>
std::enable_if_t<!std::is_same_v<COEFF,int>, Constraint<COEFF>>
operator <=(int c, const LinearSum<COEFF> &linExp) {
    return linExp >= c;
}

template<class COEFF>
std::enable_if_t<!std::is_same_v<COEFF,int>, Constraint<COEFF>>
operator >=(int c, const LinearSum<COEFF> &linExp) {
    return linExp <= c;
}


template<class COEFF>
std::ostream &operator <<(std::ostream &out, const Constraint<COEFF> &constraint) {
    out << constraint.lowerBound << " <= ";
    bool first = true;
    for(int term=0; term < constraint.coefficients.sparseSize(); ++term) {
        if(first) first = false; else out << " + \t";
        out << constraint.coefficients.values[term] << "X(" << constraint.coefficients.indices[term] << ")";
    }
    out << " <= " << constraint.upperBound;
    return  out;
}


#endif //GLPKPP_CONSTRAINT_H
