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

#ifndef ABMCMC_SPARSEFUNCTION_H
#define ABMCMC_SPARSEFUNCTION_H

#include <vector>
#include <array>
#include <functional>
#include <set>
#include "subscript_operator_traits.h"

template<typename OUT, typename IN>
class SparseFunction: public std::function<OUT(IN)> {
public:
    typedef typename subscript_operator_traits<IN>::base_type element_type;

    std::vector<int>          dependencies;


    SparseFunction(std::function<OUT(IN)> widenedFunction, std::vector<int> dependencies):
            std::function<OUT(IN)>(std::move(widenedFunction)),
            dependencies(std::move(dependencies))
    { }


    // Convenience constructors for functions with low numbers of arguments
    SparseFunction(std::function<OUT(const element_type &)> unaryFunction, int argumentIndex):
            std::function<OUT(IN)>([unaryFunction, argumentIndex](IN args) {
                return unaryFunction(args[argumentIndex]);
            }),
            dependencies({argumentIndex})
    { }

    SparseFunction(std::function<OUT(const element_type &, const element_type &)> binaryFunction, int argumentIndex1, int argumentIndex2):
            std::function<OUT(IN)>([binaryFunction, argumentIndex1, argumentIndex2](IN args) {
                return binaryFunction(args[argumentIndex1], args[argumentIndex2]);
            }),
            dependencies({argumentIndex1, argumentIndex2})
    { }


    // check whether non-dependent vars are truly non-dependent
    void sanityCheck(std::remove_const_t<std::remove_reference_t<IN>> inputVector) const {
        std::set<int> deps(dependencies.begin(), dependencies.end());
        OUT baseValue = (*this)(inputVector);
        for(int t=0; t<100; ++t) {
            for(int i=0; i< inputVector.size(); ++i) {
                if(deps.find(i) == deps.end()) {
                    inputVector[i] += Random::nextBool()?-1:1;
                }
                assert((*this)(inputVector) == baseValue);
            }
        }
    }

//    OUT widenedValue(IN X) const { return function(X).first; }
//
//    OUT exactValue(IN X) const {
//        auto cValue = function(X);
//        return cValue.second?cValue.first:-INFINITY;
//    }
};


#endif //ABMCMC_SPARSEFUNCTION_H
