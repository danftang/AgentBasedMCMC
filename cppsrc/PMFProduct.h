//
// Created by daniel on 10/08/2021.
//

#ifndef GLPKTEST_PMFPRODUCT_H
#define GLPKTEST_PMFPRODUCT_H

// represents the product of a number of PMFs
class PMFProduct {
public:
    std::vector<std::function<double(const std::vector<double> &)>> pmfs;

    int size() const { return pmfs.size(); }

    double logP(const std::vector<double> &X) const {
        double lP = 0.0;
        for(const std::function<double(const std::vector<double> &)> &pmf: pmfs) lP += pmf(X);
        return lP;
    }

    PMFProduct &operator *=(std::function<double(const std::vector<double> &)> pmf) {
        pmfs.push_back(std::move(pmf));
        return *this;
    }

    PMFProduct &operator *=(PMFProduct others) {
        for(std::function<double(const std::vector<double> &)> &pmf: others.pmfs) pmfs.push_back(std::move(pmf));
        return *this;
    }

    PMFProduct operator *(std::function<double(const std::vector<double> &)> pmf) const & {
        PMFProduct prod;
        prod.pmfs.reserve(size() + 1);
        for(const std::function<double(const std::vector<double> &)> &thispmf: pmfs) prod.pmfs.push_back(thispmf);
        prod *= pmf;
        return prod;
    }

    PMFProduct operator *(std::function<double(const std::vector<double> &)> pmf) && {
        (*this) *= pmf;
        return std::move(*this);
    }

    std::function<double(const std::vector<double> &)> PMF() const & {
        return [*this](const std::vector<double> &X) { return logP(X); };
    }

    std::function<double(const std::vector<double> &)> PMF() && {
        return [prod = std::move(*this)](const std::vector<double> &X) { return prod.logP(X); };
    }


};


#endif //GLPKTEST_PMFPRODUCT_H
