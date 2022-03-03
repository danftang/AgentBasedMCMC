//
// Created by daniel on 10/08/2021.
//

#ifndef GLPKTEST_PMFPRODUCT_H
#define GLPKTEST_PMFPRODUCT_H

// represents the product of a number of PMFs
template<typename DOMAIN>
class PMFProduct {
public:
    std::vector<std::function<double(const DOMAIN &)>> pmfs;

    int size() const { return pmfs.size(); }

    double operator ()(const DOMAIN &X) const { return logP(X); }
    double logP(const DOMAIN &X) const {
        double lP = 0.0;
        for(const std::function<double(const DOMAIN &)> &pmf: pmfs) lP += pmf(X);
        return lP;
    }

    PMFProduct &operator *=(std::function<double(const DOMAIN &)> pmf) {
        pmfs.push_back(std::move(pmf));
        return *this;
    }

    PMFProduct &operator *=(PMFProduct others) {
        for(std::function<double(const DOMAIN &)> &pmf: others.pmfs) pmfs.push_back(std::move(pmf));
        return *this;
    }

    PMFProduct operator *(std::function<double(const DOMAIN &)> pmf) const & {
        PMFProduct prod;
        prod.pmfs.reserve(size() + 1);
        for(const std::function<double(const DOMAIN &)> &thispmf: pmfs) prod.pmfs.push_back(thispmf);
        prod *= pmf;
        return prod;
    }

    PMFProduct operator *(std::function<double(const DOMAIN &)> pmf) && {
        (*this) *= pmf;
        return std::move(*this);
    }

//    std::function<double(const DOMAIN &)> PMF() const & {
//        return [*this](const DOMAIN &X) { return logP(X); };
//    }
//
//    std::function<double(const std::vector<double> &)> PMF() && {
//        return [prod = std::move(*this)](const std::vector<double> &X) { return prod.logP(X); };
//    }


};


#endif //GLPKTEST_PMFPRODUCT_H
