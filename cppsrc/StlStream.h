//
// Created by daniel on 07/05/2021.
//

#include <ostream>
#include <vector>
#include <iterator>

#ifndef GLPKTEST_STLSTREAM_H
#define GLPKTEST_STLSTREAM_H

template<typename T>
std::ostream &operator <<(std::ostream &out, const std::vector<T> &vec) {
    out << "{";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b}";
    return out;
}

#endif //GLPKTEST_STLSTREAM_H
