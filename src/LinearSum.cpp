//
// Created by daniel on 27/04/2021.
//
#include <iomanip>
#include <iostream>

#include "LinearSum.h"
#include "Constraint.h"

    Constraint operator ==(const LinearSum &linExp, double c) {
        return Constraint(c, linExp.toSparseVec(), c);
    }

    Constraint operator <=(const LinearSum &linExp, double c) {
        return Constraint(-std::numeric_limits<double>::infinity(), linExp.toSparseVec(), c);
    }

    Constraint operator >=(const LinearSum &linExp, double c) {;
        return Constraint(c, linExp.toSparseVec(), std::numeric_limits<double>::infinity());
    }

    Constraint operator ==(double c, const LinearSum &linExp) {
        return linExp == c;
    }

    Constraint operator <=(double c, const LinearSum &linExp) {
        return linExp >= c;
    }

    Constraint operator >=(double c, const LinearSum &linExp) {
        return linExp <= c;
    }

    Constraint operator==(const LinearSum::Monomial &monomial, double c) {
        Constraint constraint(c,c);
        constraint.coefficients.insert(monomial.variableId, monomial.multiplier);
        return constraint;
    }

    Constraint operator<=(const LinearSum::Monomial &monomial, double c) {
        Constraint constraint(-std::numeric_limits<double>::infinity(),c);
        constraint.coefficients.insert(monomial.variableId, monomial.multiplier);
        return constraint;
    }

    Constraint operator>=(const LinearSum::Monomial &monomial, double c) {
        Constraint constraint(c,std::numeric_limits<double>::infinity());
        constraint.coefficients.insert(monomial.variableId, monomial.multiplier);
        return constraint;
    }

    Constraint operator==(double c, const LinearSum::Monomial &monomial) {
        return monomial == c;
    }

    Constraint operator<=(double c, const LinearSum::Monomial &monomial) {
        return monomial >= c;
    }

    Constraint operator>=(double c, const LinearSum::Monomial &monomial) {
        return monomial <= c;
    }

    std::ostream &operator <<(std::ostream &out, const LinearSum &sum) {
        bool first = true;
        for (auto [varId, coeff] : sum.coefficients) {
            if(first) first = false; else out << " +\t";
            out << std::setw(8) << coeff << "X" << varId;
        }
        return out;
    }


