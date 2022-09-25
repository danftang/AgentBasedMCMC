// Created by daniel on 02/08/2021.
//

#ifndef GLPKTEST_ABMCONSTRAINTS_H
#define GLPKTEST_ABMCONSTRAINTS_H

class ABM {
public:
//    typedef int occupation_type;    // data type for act- and state-occupation numbers
//    typedef int coefficient_type;   // data type for constraints
    static thread_local double kappa;
};

inline thread_local double ABM::kappa = 4.0;


#endif //GLPKTEST_ABMCONSTRAINTS_H
