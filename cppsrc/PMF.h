//
// Created by daniel on 26/07/2021.
//

#ifndef GLPKTEST_PMF_H
#define GLPKTEST_PMF_H


class PMF {
public:
    virtual double              logProb(const std::vector<double> &X) = 0;
    virtual std::vector<double> nextSample() = 0;

    virtual ~PMF();
};


#endif //GLPKTEST_PMF_H
