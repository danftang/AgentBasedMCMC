//
// Created by daniel on 01/09/22.
//

#ifndef ABMCMC_EQUALITYCONSTRAINT_H
#define ABMCMC_EQUALITYCONSTRAINT_H

#include "include/SparseVec.h"
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

template<class CONST, class COEFF, class RESULT = decltype(std::declval<COEFF>() - std::declval<CONST>())>
EqualityConstraint<RESULT> operator ==(const LinearSum<COEFF> &linExp, CONST c) {
    EqualityConstraint<RESULT> result;
    result.constant = c;
    result.coefficients.reserve(linExp.coefficients.size());
    for(const auto &entry: linExp.coefficients) result.coefficients.insert(entry.first, entry.second);
    return result;
}

template<class COEFF, class CONST, class RESULT = decltype(std::declval<COEFF>() - std::declval<CONST>())>
EqualityConstraint<COEFF> operator ==(CONST c, const LinearSum<COEFF> &linExp) {
    return operator ==<CONST,COEFF,RESULT>(linExp,c);
}

template<class COEFF, typename = std::enable_if_t<std::is_arithmetic_v<COEFF>>>
EqualityConstraint<COEFF> operator ==(const X &symbolicVar, const COEFF &constant) {
    return EqualityConstraint<COEFF>({{symbolicVar.id,1}}, constant);
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
