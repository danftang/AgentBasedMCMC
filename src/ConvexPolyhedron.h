//
// Created by daniel on 30/07/2021.
//

#ifndef GLPKTEST_CONVEXPOLYHEDRON_H
#define GLPKTEST_CONVEXPOLYHEDRON_H


#include <vector>
#include "StlStream.h"
//#include "agents/PredPreyAgent.h"
#include "Constraint.h"

class ConvexPolyhedron: public std::vector<Constraint> {
public:
    ConvexPolyhedron(): std::vector<Constraint>() {}
    explicit ConvexPolyhedron(const std::vector<Constraint> &constraints): std::vector<Constraint>(constraints) {}
    explicit ConvexPolyhedron(std::vector<Constraint> &&constraints): std::vector<Constraint>(constraints) {}
    ConvexPolyhedron(const std::initializer_list<Constraint> &constraints): std::vector<Constraint>(constraints) {}

    bool isValidSolution(const std::vector<double> &X) const {
        for(const Constraint &constraint: *this) {
//            std::cout << "Checking if " << X << " satisfies " << constraint << std::endl;
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
        for(const Constraint &constraint: *this) {
            dim = std::max(dim,constraint.coefficients.maxNonZeroIndex());
        }
        return dim;
    }

    ConvexPolyhedron &operator +=(const std::vector<Constraint> &other) {
        reserve(size()+other.size());
        for(const Constraint &constraint: other) {
            push_back(constraint);
        }
        return *this;
    }

    ConvexPolyhedron &operator +=(std::vector<Constraint> &&other) {
        reserve(size()+other.size());
        for(Constraint &constraint: other) {
            push_back(std::move(constraint));
        }
        return *this;
    }

    ConvexPolyhedron operator +(const std::vector<Constraint> &other) && {
        (*this) += other;
        return std::move(*this);
    }

    ConvexPolyhedron operator +(std::vector<Constraint> &&other) && {
        (*this) += std::move(other);
        return std::move(*this);
    }

    ConvexPolyhedron operator +(std::vector<Constraint> &&other) const & {
        ConvexPolyhedron result(std::move(other));
        result += *this;
        return result;
    }

    ConvexPolyhedron operator +(const std::vector<Constraint> &other) const & {
        ConvexPolyhedron result(*this);
        result += other;
        return result;
    }

    friend std::ostream &operator <<(std::ostream &out, const ConvexPolyhedron &convexPolyhedron) {
        for(const Constraint &constraint: convexPolyhedron) {
            out << constraint << std::endl;
        }
        return out;
    }

};


#endif //GLPKTEST_CONVEXPOLYHEDRON_H
