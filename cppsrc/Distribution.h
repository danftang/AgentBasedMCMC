//
// Created by daniel on 09/08/2021.
//

#ifndef GLPKTEST_DISTRIBUTION_H
#define GLPKTEST_DISTRIBUTION_H

// Before we know whether we want to sample from or take the log prob of this
// distribution, then we can use the more abstract Distribution object.
class Distribution {
public:
    typedef std::function<std::vector<double>()>              Sampler;
//    typedef std::function<double(const std::vector<double> &)>  LogPMF;

    virtual ConvexPMF   PMF() const =0;
    virtual Sampler     sampler() const =0;
    virtual int         dimension() const =0;
};


#endif //GLPKTEST_DISTRIBUTION_H
