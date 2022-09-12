//
// Created by daniel on 01/09/22.
//

#ifndef ABMCMC_EQUALITYCONSTRAINT_H
#define ABMCMC_EQUALITYCONSTRAINT_H

#include "SparseVec.h"
#include "LinearSum.h"

template<class COEFF>
class EqualityConstraint {
public:

    SparseVec<COEFF>    coefficients;
    COEFF               constant;

    EqualityConstraint(): EqualityConstraint(0) {} // default to \emptyset = 0
    explicit EqualityConstraint(COEFF constant): constant(constant) {};
    EqualityConstraint(SparseVec<COEFF> coefficients, COEFF constant): coefficients(std::move(coefficients)), constant(constant) {};

    template<typename T>
    bool isValidSolution(const T &X) const {
        return coefficients * X == constant;
    }

    int maxCoefficientIndex() const {
        if(coefficients.indices.size() == 0) return 0;
        return *std::max_element(coefficients.indices.begin(), coefficients.indices.end());
    }

};

template<class COEFF>
EqualityConstraint<COEFF> operator ==(const LinearSum<COEFF> &linExp, COEFF c) {
    return EqualityConstraint<COEFF>(linExp.toSparseVec(), c);
}

template<class COEFF>
EqualityConstraint<COEFF> operator ==(double c, const LinearSum<COEFF> &linExp) {
    return linExp == c;
}

template<typename COEFF>
std::ostream &operator <<(std::ostream &out, const EqualityConstraint<COEFF> &constraint) {
    for(int i=0; i < constraint.coefficients.sparseSize(); ++i) {
        out << constraint.coefficients.values[i] << "X" << constraint.coefficients.indices[i];
        if(i != constraint.coefficients.sparseSize()-1) out << " + ";
    }
    out << " == " << constraint.constant;
    return out;
}



#endif //ABMCMC_EQUALITYCONSTRAINT_H
