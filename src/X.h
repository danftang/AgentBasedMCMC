// Represents a variable X_i for use with Constraint, LienarSum and ConvexPolyhedron
//
// Created by daniel on 26/04/2021.
//

#ifndef ABMCMC_X_H
#define ABMCMC_X_H

class X {
public:
    int id;
    explicit X(int id): id(id) { }

//    operator int() const { return id; }
};

#endif

