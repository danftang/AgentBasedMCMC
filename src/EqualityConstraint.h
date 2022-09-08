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

//    EqualityConstraint(): EqualityConstraint(0) {}
    explicit EqualityConstraint(COEFF constant): constant(constant) {};
    EqualityConstraint(SparseVec<COEFF> coefficients, COEFF constant): coefficients(std::move(coefficients)), constant(constant) {};

    bool isValidSolution(const std::vector<COEFF> &X) const {
        return coefficients * X == constant;
    }

    friend std::ostream &operator <<(std::ostream &out, const EqualityConstraint<COEFF> &constraint) {
        for(int i=0; i < constraint.coefficients.sparseSize(); ++i) {
            out << constraint.coefficients.values[i] << "X" << constraint.coefficients.indices[i];
            if(i != constraint.coefficients.sparseSize()-1) out << " + ";
        }
        out << " == " << constraint.constant;
        return out;
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



#endif //ABMCMC_EQUALITYCONSTRAINT_H
