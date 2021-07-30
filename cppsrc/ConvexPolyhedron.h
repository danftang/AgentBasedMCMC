//
// Created by daniel on 30/07/2021.
//

#ifndef GLPKTEST_CONVEXPOLYHEDRON_H
#define GLPKTEST_CONVEXPOLYHEDRON_H


#include <vector>
#include "glpkpp.h"

class ConvexPolyhedron: public std::vector<glp::Constraint> {
public:

    ConvexPolyhedron &operator +=(const ConvexPolyhedron &other) {
        reserve(size()+other.size());
        for(const glp::Constraint &constraint: other) {
            push_back(constraint);
        }
        return *this;
    }

    ConvexPolyhedron &operator +=(ConvexPolyhedron &&other) {
        reserve(size()+other.size());
        for(glp::Constraint &constraint: other) {
            push_back(std::move(constraint));
        }
        return *this;
    }


    ConvexPolyhedron operator +(const ConvexPolyhedron &other) && {
        std::cout << "Doung lhs move\n";
        (*this) += other;
        return std::move(*this);
    }

    ConvexPolyhedron operator +(ConvexPolyhedron &&other) && {
        std::cout << "Doung double move\n";
        (*this) += std::move(other);
        return std::move(*this);
    }

    ConvexPolyhedron operator +(ConvexPolyhedron &&other) const & {
        std::cout << "Doung rhs move\n";
        other += *this;
        return std::move(other);
    }

    ConvexPolyhedron operator +(const ConvexPolyhedron &other) const & {
        std::cout << "Doung plain addition\n";
        ConvexPolyhedron result(*this);
        result += other;
        return result;
    }

    friend std::ostream &operator <<(std::ostream &out, const ConvexPolyhedron &convexPolyhedron) {
        for(const glp::Constraint &constraint: convexPolyhedron) {
            out << constraint << std::endl;
        }
        return out;
    }

};


#endif //GLPKTEST_CONVEXPOLYHEDRON_H
