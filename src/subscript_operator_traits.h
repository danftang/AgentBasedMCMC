//
// Created by daniel on 09/09/22.
//

#ifndef ABMCMC_SUBSCRIPT_OPERATOR_TRAITS_H
#define ABMCMC_SUBSCRIPT_OPERATOR_TRAITS_H

#include <type_traits>
#include <functional>

// If T is a type that has a subscript operator that takes an integer T[i]
// then subsectipt_operator_type<T>::return_type is the return type of T[i]
// and  subsectipt_operator_type<T>::base_type is return_type with any references and cv qualification removed


template<typename T, typename = decltype(std::declval<T>()[0])>
class subscript_operator_traits {
public:
    typedef decltype(std::declval<T>()[0]) return_type;
    typedef std::remove_cv_t<std::remove_reference_t<return_type>> base_type;
};


#endif //ABMCMC_SUBSCRIPT_OPERATOR_TRAITS_H
