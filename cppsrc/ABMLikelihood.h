//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_ABMLIKELIHOOD_H
#define GLPKTEST_ABMLIKELIHOOD_H

template<typename AGENT>
class ABMLikelihood {
    std::map<State<AGENT>,std::function<double(int)>> likelihoods;
    ConvexPolyhedron constraints;


};


#endif //GLPKTEST_ABMLIKELIHOOD_H
