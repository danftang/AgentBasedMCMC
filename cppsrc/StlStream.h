//
// Created by daniel on 07/05/2021.
//

#include <ostream>
#include <vector>
#include <iterator>
#include <map>

#ifndef GLPKTEST_STLSTREAM_H
#define GLPKTEST_STLSTREAM_H

template<typename T>
std::ostream &operator <<(std::ostream &out, const std::vector<T> &vec) {
    out << "{";
//    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(out, ", "));
    for(int i=0; i<vec.size(); ++i) out << vec[i] << ", ";
    out << "\b\b}";
    return out;
}

template<typename T>
std::ostream &operator <<(std::ostream &out, const std::valarray<T> &vec) {
    out << "{";
    for(int i=0; i<vec.size(); ++i) out << vec[i] << ", ";
    out << "\b\b}";
    return out;
}


template<typename KEY, typename VALUE>
std::ostream &operator <<(std::ostream &out, const std::map<KEY,VALUE> &map) {
    out << "{";
    for(const auto &[key,value]: map) {
        out << key << " -> " << value << ", ";
    }
    out << "\b\b}";
    return out;
}

template<typename KEY, typename VALUE>
std::ostream &operator <<(std::ostream &out, const std::multimap<KEY,VALUE> &map) {
    out << "{";
    for(const auto &[key,value]: map) {
        out << key << " -> " << value << ", ";
    }
    out << "\b\b}";
    return out;
}


template<typename T1, typename T2>
std::ostream &operator <<(std::ostream &out, const std::pair<T1,T2> &pair) {
    out << "(" << pair.first << ", " << pair.second << ")";
    return out;
}

#endif //GLPKTEST_STLSTREAM_H
