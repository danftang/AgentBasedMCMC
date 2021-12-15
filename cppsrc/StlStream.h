//
// Created by daniel on 07/05/2021.
//

#include <ostream>
#include <vector>
#include <list>
#include <iterator>
#include <map>
#include <chrono>

#ifndef GLPKTEST_STLSTREAM_H
#define GLPKTEST_STLSTREAM_H

//template <typename T, typename=void>
//static constexpr bool hasBeginEnd = false;


//template <typename T>
//static constexpr bool hasBeginEnd<T, std::void_t<
//        decltype(std::declval<T>().begin()),
//        decltype(std::declval<T>().end())
//>> = true;

template<typename T> std::ostream & operator <<(std::ostream &out, const std::list<T> &container);

template<typename T>
std::ostream & operator <<(std::ostream &out, const std::vector<T> &container) {
    out << "{";
    typename std::vector<T>::const_iterator it = container.begin();
    if(it != container.end()) {
        out << *it;
        while(++it != container.end()) out << ", " << *it;
    }
    out << "}";
    return out;
}

template<typename T>
std::ostream & operator <<(std::ostream &out, const std::list<T> &container) {
    out << "{";
    auto it = container.begin();
    if(it != container.end()) {
        out << *it;
        while(++it != container.end()) out << ", " << *it;
    }
    out << "}";
    return out;
}




//template<typename T>
//std::ostream &operator <<(std::ostream &out, const std::vector<T> &vec) {
//    out << "{";
//    for(int i=0; i<vec.size()-1; ++i) out << vec[i] << ", ";
//    if(vec.size() > 0) out << vec[vec.size()-1];
//    out << "}";
//    return out;
//}


template<typename T>
std::ostream &operator <<(std::ostream &out, const std::valarray<T> &vec) {
    out << "{";
    for(int i=0; i<vec.size()-1; ++i) out << vec[i] << ", ";
    if(vec.size() > 0) out << vec[vec.size()-1];
    out << "}";
    return out;
}


template<typename KEY, typename VALUE>
std::ostream &operator <<(std::ostream &out, const std::map<KEY,VALUE> &map) {
    out << "{";
    for(const auto &[key,value]: map) {
        out << key << " -> " << value << " ";
    }
    out << "}";
    return out;
}

template<typename KEY, typename VALUE>
std::ostream &operator <<(std::ostream &out, const std::multimap<KEY,VALUE> &map) {
    out << "{";
    for(const auto &[key,value]: map) {
        out << key << " -> " << value << " ";
    }
    out << "}";
    return out;
}


template<typename T1, typename T2>
std::ostream &operator <<(std::ostream &out, const std::pair<T1,T2> &pair) {
    out << "(" << pair.first << ", " << pair.second << ")";
    return out;
}

template<typename I, intmax_t UNITS>
std::ostream &operator <<(std::ostream &out, const std::chrono::duration<I,std::ratio<1,UNITS>> duration) {
    long double count = duration.count();
    intmax_t units = UNITS;
    while(count >= 1000.0 && units >= 1000) {
        count /= 1000.0;
        units /= 1000;
    }
    out << count;
    switch(units) {
        case 1: out << "s"; break;
        case 1000: out << "ms"; break;
        case 1000000: out << "μs"; break;
        case 1000000000: out << "ns"; break;
    }
    return out;
}

#endif //GLPKTEST_STLSTREAM_H
