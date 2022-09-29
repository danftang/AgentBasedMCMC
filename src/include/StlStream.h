//
// Created by daniel on 07/05/2021.
//

#include <ostream>
#include <vector>
#include <list>
#include <iterator>
#include <map>
#include <set>
#include <chrono>
#include <deque>
#include <forward_list>
#include <unordered_set>

#ifndef GLPKTEST_STLSTREAM_H
#define GLPKTEST_STLSTREAM_H

// This class allows us to implement a version of operator << that requires an implicit conversion
// of its arguments. This ensures that a class that has begin() and end() operators has a default
// stream operator but any direct implementation will override the default.
//struct ostream_implicit_conversion {
//    std::ostream &out;
//    ostream_implicit_conversion(std::ostream &out) : out(out) {}
//};
//
//template<typename T, typename = decltype(std::declval<T>().begin()), typename = decltype(std::declval<T>().end())>
//std::ostream & operator <<(ostream_implicit_conversion out, const T &container) {
//    out.out << "{";
//    auto it = container.begin();
//    if(it != container.end()) {
//        out.out << *it;
//        while(++it != container.end()) out.out << ", " << *it;
//    }
//    out.out << "}";
//    return out.out;
//}

template<class T>
std::ostream & default_container_printer(std::ostream &out, const T &container) {
    out << "{";
    auto it = container.begin();
    if(it != container.end()) {
        out << *it;
        while(++it != container.end()) out << ", " << *it;
    }
    out << "}";
    return out;
}

template<class T, int N>
std::ostream & operator <<(std::ostream &out, const std::array<T,N> &container) { return default_container_printer(out, container); }
template<class T>
std::ostream & operator <<(std::ostream &out, const std::vector<T> &container) { return default_container_printer(out, container); }
template<class T>
std::ostream & operator <<(std::ostream &out, const std::deque<T> &container) { return default_container_printer(out, container); }
template<class T>
std::ostream & operator <<(std::ostream &out, const std::forward_list<T> &container) { return default_container_printer(out, container); }
template<class T>
std::ostream & operator <<(std::ostream &out, const std::list<T> &container) { return default_container_printer(out, container); }
template<class T>
std::ostream & operator <<(std::ostream &out, const std::set<T> &container) { return default_container_printer(out, container); }
template<class T>
std::ostream & operator <<(std::ostream &out, const std::unordered_set<T> &container) { return default_container_printer(out, container); }



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



template<typename I, intmax_t UNITS>
std::ostream &operator <<(std::ostream &out, const std::chrono::duration<I,std::ratio<1,UNITS>> &duration) {
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

template<typename T1, typename T2>
std::ostream &operator <<(std::ostream &out, const std::pair<T1,T2> &pair) {
    out << "(" << pair.first << ", " << pair.second << ")";
    return out;
}


#endif //GLPKTEST_STLSTREAM_H
