//
// Created by daniel on 09/09/22.
//

#ifndef ABMCMC_EXTRATRAITS_H
#define ABMCMC_EXTRATRAITS_H

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


template<class T, typename Head, typename... Tail>
class tuple_index {
    static constexpr int value = tuple_index<T,Tail...>::value + 1;
};

template<class T, typename... Tail>
class tuple_index<T,T,Tail...> {
    static constexpr int value = 0;
};

template<class T, class H>
class tuple_index<T,H> {
    static constexpr int value = -1;
};

template<class T>
class tuple_index<T,T> {
    static constexpr int value = 0;
};


template<class T>
class subset_index: public tuple_index<T, bool, char, short, int, long, float, double> {};

template<> class subset_index<int>: public tuple_index<int, bool, char, unsigned char, short, unsigned short, int> {};

template<class T>
inline static constexpr int subset_index_v = subset_index<T>::value;


template<class A, class B, bool V = subset_index_v<A> != -1 && subset_index_v<B> != -1 &&  subset_index_v<A> <= subset_index_v<B>>
class is_subset_of {
    static constexpr bool value = V;
};






#endif //ABMCMC_EXTRATRAITS_H
