// Represents a pair of functions, Fe and Fw, from an ARGS to an OUT
// where ARGS should have a subscript operator.
// They are sparse in the sense that their value depends only on
// a small fraction of elements of the input vector, and are
// independent of the rest. the `depenencies' vector defines which
// elements of X the functions are dependent on.
//
// Fw is a "widened" version of Fe in the sense that
// if Fe is not -infinity then Fw equals Fe
// however, if Fe is -infinity, then Fw could be finite.
// Both functions are defined by `widenedFunction' which
// returns a pair containing the value of Fw and a boolean
// which is true if Fe is not -infinity.
//
// TODO: Should this be generalised to importance sampling?
//       where Fw is any function that is never smaller than
//       Fe.
//
// Created by daniel on 01/09/22.
//

#ifndef ABMCMC_SPARSEWIDENEDFUNCTION_H
#define ABMCMC_SPARSEWIDENEDFUNCTION_H

#include <vector>
#include <array>
#include <functional>

template<typename OUT, typename ARGS, typename ARG_ELEMENT = decltype(std::declval<ARGS>()[0])>
class SparseWidenedFunction {
public:
    std::function<std::pair<OUT,bool>(ARGS)>    widenedFunction;
    std::vector<int>                            dependencies;

    SparseWidenedFunction(std::function<std::pair<OUT,bool>(ARGS)> widenedFunction, std::initializer_list<int> dependencies):
        widenedFunction(std::move(widenedFunction)),
        dependencies(dependencies)
    {
    }

    SparseWidenedFunction(std::function<std::pair<OUT,bool>(ARGS)> widenedFunction, std::vector<int> dependencies):
            widenedFunction(std::move(widenedFunction)),
            dependencies(std::move(dependencies))
    { }


    // Convenience constructors for functions with low numbers of arguments
    SparseWidenedFunction(std::function<std::pair<OUT,bool>(const ARG_ELEMENT &)> widenedUnaryFunction, int argumentIndex):
            widenedFunction([widenedUnaryFunction, argumentIndex](ARGS args) {
                return widenedUnaryFunction(args[argumentIndex]);
            }),
            dependencies({argumentIndex})
    { }

    SparseWidenedFunction(std::function<std::pair<OUT,bool>(const ARG_ELEMENT &, const ARG_ELEMENT &)> widenedUnaryFunction, int argumentIndex1, int argumentIndex2):
            widenedFunction([widenedUnaryFunction, argumentIndex1, argumentIndex2](ARGS args) {
                return widenedUnaryFunction(args[argumentIndex1], args[argumentIndex2]);
            }),
            dependencies({argumentIndex1, argumentIndex2})
    { }


    OUT widenedValue(ARGS X) const { return widenedFunction(X).first; }

    OUT exactValue(ARGS X) const {
        auto cValue = widenedFunction(X);
        return cValue.second?cValue.first:-INFINITY;
    }

};


#endif //ABMCMC_SPARSEWIDENEDFUNCTION_H
