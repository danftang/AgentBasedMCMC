//
// Created by daniel on 23/03/2022.
//

#ifndef ABMCMC_PERTURBABLEFUNCTION_H
#define ABMCMC_PERTURBABLEFUNCTION_H

#include "include/SparseVec.h"

template<class D, class RANGE>
class PerturbableFunction {
public:
    virtual ~PerturbableFunction() {};

    virtual void perturb(const std::vector<D> &stateBeforePerturbation, const std::vector<int> &perturbedIndices)=0;
    virtual void undo()=0; // returns an error if last perturbation was not undoable
    virtual RANGE getLogValue(const std::vector<D> &currentState)=0;
    virtual void setState(const std::vector<D> &currentState)=0; // without assuming any elements are unchanged
};


#endif //ABMCMC_PERTURBABLEFUNCTION_H
