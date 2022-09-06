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

};


#endif //ABMCMC_EQUALITYCONSTRAINTS_H
