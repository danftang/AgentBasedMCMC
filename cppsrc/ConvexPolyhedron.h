//
// Created by daniel on 30/07/2021.
//

#ifndef GLPKTEST_CONVEXPOLYHEDRON_H
#define GLPKTEST_CONVEXPOLYHEDRON_H


#include <vector>
#include "glpkpp.h"

class ConvexPolyhedron: public std::vector<glp::Constraint> {
public:
    ConvexPolyhedron(): std::vector<glp::Constraint>() {}
    explicit ConvexPolyhedron(const std::vector<glp::Constraint> &constraints): std::vector<glp::Constraint>(constraints) {}
    explicit ConvexPolyhedron(std::vector<glp::Constraint> &&constraints): std::vector<glp::Constraint>(constraints) {}
    ConvexPolyhedron(const std::initializer_list<glp::Constraint> &constraints): std::vector<glp::Constraint>(constraints) {}

    bool isValidSolution(const std::vector<double> &X) const {
        for(const glp::Constraint &constraint: *this) {
            if(!constraint.isValidSolution(X)) return false;
        }
        return true;

    }

    glp::Problem toLPProblem() const {
        std::cout << "Constraints are:\n" << *this << std::endl;
        glp::Problem lp(*this);
        lp.advBasis();
        lp.warmUp();
        std::cout << "Problem is:\n" << lp << std::endl;
        return lp;
    }


    ConvexPolyhedron &operator +=(const std::vector<glp::Constraint> &other) {
        reserve(size()+other.size());
        for(const glp::Constraint &constraint: other) {
            push_back(constraint);
        }
        return *this;
    }

    ConvexPolyhedron &operator +=(std::vector<glp::Constraint> &&other) {
        reserve(size()+other.size());
        for(glp::Constraint &constraint: other) {
            push_back(std::move(constraint));
        }
        return *this;
    }

    ConvexPolyhedron operator +(const std::vector<glp::Constraint> &other) && {
        (*this) += other;
        return std::move(*this);
    }

    ConvexPolyhedron operator +(std::vector<glp::Constraint> &&other) && {
        (*this) += std::move(other);
        return std::move(*this);
    }

    ConvexPolyhedron operator +(std::vector<glp::Constraint> &&other) const & {
        ConvexPolyhedron result(std::move(other));
        result += *this;
        return result;
    }

    ConvexPolyhedron operator +(const std::vector<glp::Constraint> &other) const & {
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
