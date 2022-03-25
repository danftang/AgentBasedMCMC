//
// Created by daniel on 30/07/2021.
//

#ifndef GLPKTEST_CONVEXPOLYHEDRON_H
#define GLPKTEST_CONVEXPOLYHEDRON_H


#include <vector>
#include <iostream>
#include "StlStream.h"
#include "Constraint.h"

template<class T>
class ConvexPolyhedron: public std::vector<Constraint<T>> {
public:
    ConvexPolyhedron(): std::vector<Constraint<T>>() {}
    ConvexPolyhedron(const std::initializer_list<Constraint<T>> &constraints): std::vector<Constraint<T>>(constraints) {}

    bool isValidSolution(const std::vector<T> &X) const {
        for(const Constraint<T> &constraint: *this) {
            if(!constraint.isValidSolution(X)) {
//                std::cout << X << " does not satisfy " << constraint << std::endl;
                return false;
            }
        }
        return true;
    }


    // the id of the highest variable in this polyhedron
    int dimension() {
        int dim = 0;
        for(const Constraint<T> &constraint: *this) {
            dim = std::max(dim,constraint.coefficients.maxNonZeroIndex());
        }
        return dim;
    }

    ConvexPolyhedron<T> &operator +=(const std::vector<Constraint<T>> &other) {
        this->reserve(this->size()+other.size());
        for(const Constraint<T> &constraint: other) {
            this->push_back(constraint);
        }
        return *this;
    }

    ConvexPolyhedron<T> &operator +=(std::vector<Constraint<T>> &&other) {
        this->reserve(this->size()+other.size());
        for(Constraint<T> &constraint: other) {
            this->push_back(std::move(constraint));
        }
        return *this;
    }

    ConvexPolyhedron<T> operator +(const std::vector<Constraint<T>> &other) && {
        (*this) += other;
        return std::move(*this);
    }

    ConvexPolyhedron<T> operator +(std::vector<Constraint<T>> &&other) && {
        (*this) += std::move(other);
        return std::move(*this);
    }

    ConvexPolyhedron<T> operator +(std::vector<Constraint<T>> &&other) const & {
        ConvexPolyhedron result(std::move(other));
        result += *this;
        return result;
    }

    ConvexPolyhedron<T> operator +(const std::vector<Constraint<T>> &other) const & {
        ConvexPolyhedron result(*this);
        result += other;
        return result;
    }

    friend std::ostream &operator <<(std::ostream &out, const ConvexPolyhedron<T> &convexPolyhedron) {
        for(const Constraint<T> &constraint: convexPolyhedron) {
            out << constraint << std::endl;
        }
        return out;
    }

};


#endif //GLPKTEST_CONVEXPOLYHEDRON_H
