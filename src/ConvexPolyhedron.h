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
    explicit ConvexPolyhedron(const std::vector<Constraint<T>> &constraints): std::vector<Constraint<T>>(constraints) {}
    explicit ConvexPolyhedron(std::vector<Constraint<T>> &&constraints): std::vector<Constraint<T>>(constraints) {}
    ConvexPolyhedron(const std::initializer_list<Constraint<T>> &constraints): std::vector<Constraint<T>>(constraints) {}

    bool isValidSolution(const std::vector<T> &X) const {
        for(const Constraint<T> &constraint: *this) {
            std::cout << "Checking if " << X << " satisfies " << constraint << std::endl;
            if(!constraint.isValidSolution(X)) return false;
        }
        return true;
    }


//    glp::Problem toLPProblem(const std::vector<double> &objective = std::vector<double>()) const {
////        std::cout << "Constraints are:\n" << *this << std::endl;
//        glp::Problem lp(*this);
//        lp.advBasis();
////        predPreyBasis(lp);
//        lp.warmUp();
//        lp.setObjective(objective);
////        std::cout << "Problem is:\n" << lp << std::endl;
//        return lp;
//    }


//    static void predPreyBasis(glp::Problem &lp) {
//        // remove all fixed vars from basis
//        int nFixed = 0;
//        for(int i = 1; i <= lp.nConstraints(); i++) {
//            if(lp.getRowType(i) == GLP_FX) {
//                ++nFixed;
//                lp.setRowStat(i, GLP_NL);
//            }
//
//        }
//        // add death as basic
//        int eventId = lp.nVars() - PredPreyAgentBase::actDomainSize() + PredPreyAgentBase::ActNames::DIE + 1;
//        while(nFixed > 0) {
//            lp.setColStat(eventId, GLP_BS);
//            --nFixed;
//            eventId -= PredPreyAgentBase::actDomainSize();
//        }
//    }



    // the id of the highest variable in this polyhedron
    int dimension() {
        int dim = 0;
        for(const Constraint<T> &constraint: *this) {
            dim = std::max(dim,constraint.coefficients.maxNonZeroIndex());
        }
        return dim;
    }

    ConvexPolyhedron<T> &operator +=(const std::vector<Constraint<T>> &other) {
        reserve(this->size()+other.size());
        for(const Constraint<T> &constraint: other) {
            push_back(constraint);
        }
        return *this;
    }

    ConvexPolyhedron<T> &operator +=(std::vector<Constraint<T>> &&other) {
        reserve(this->size()+other.size());
        for(Constraint<T> &constraint: other) {
            push_back(std::move(constraint));
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
