//
// Created by daniel on 09/08/2021.
//

#ifndef GLPKTEST_DISTRIBUTION_H
#define GLPKTEST_DISTRIBUTION_H


class Distribution {
public:
    typedef std::function<std::vector<double>()>              Sampler;
//    typedef std::function<double(const std::vector<double> &)>  PMF;

    virtual ConvexPMF   PMF() const =0;
    virtual Sampler     sampler() const =0;
    virtual int         dimension() const =0;
};


#endif //GLPKTEST_DISTRIBUTION_H
