//
// Created by daniel on 01/09/22.
//

#ifndef ABMCMC_EQUALITYCONSTRAINTS_H
#define ABMCMC_EQUALITYCONSTRAINTS_H

#include <vector>
#include "EqualityConstraint.h"

template<class T>
class EqualityConstraints: public std::vector<EqualityConstraint<T>> {
public:
    EqualityConstraints<T> &operator +=(const EqualityConstraints<T> &other) {
        this->insert(this->end(), other.begin(), other.end());
        return *this;
    }
};


#endif //ABMCMC_EQUALITYCONSTRAINTS_H
