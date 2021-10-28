//
// Created by daniel on 20/10/2021.
//

#ifndef GLPKTEST_MEANANDVARIANCE_H
#define GLPKTEST_MEANANDVARIANCE_H


#include <valarray>

// DATA can be any object with .size() and operator[]
//template<typename DATA, typename = std::void_t<decltype(std::declval<DATA>()[0]), decltype(std::declval<DATA>().size())>>
class MeanAndVariance {
public:
    int nSamples;
    std::valarray<double> sum;
    std::valarray<double> sumOfSquares;

    MeanAndVariance(): nSamples(0) { }

//    MeanAndVariance(int dimension): nSamples(0), sum(0.0, dimension), sumOfSquares(0.0, dimension) { }

    bool operator()(const std::vector<double> &sample) {
        if(nSamples == 0) initArrays(sample.size());
        if(sample.size() != sum.size()) throw("Dimension of sample does not match this logger.");
        for(int i=0; i<sample.size(); ++i) {
            double xi = sample[i];
            sum[i] += xi;
            sumOfSquares[i] += xi*xi;
        }
        ++nSamples;
        return true;
    }

    auto consumer() {
        return [this](const std::vector<double> &sample) { return (*this)(sample); };
    }

    std::valarray<double> mean() const {
        return sum * (1.0/nSamples);
    }


    std::valarray<double> sampleVariance() const {
        return (sumOfSquares - sum*sum * (1.0/nSamples))*(1.0/(nSamples-1));
    }

    friend std::ostream &operator <<(std::ostream &out, const MeanAndVariance &mv) {
        out << std::endl;
        out << "means = " << mv.mean() << std::endl;
        out << "variances = " << mv.sampleVariance() << std::endl;
        return out;
    }

private:

    void initArrays(int size) {
        sum.resize(size, 0.0);
        sumOfSquares.resize(size, 0.0);
    }
};


#endif //GLPKTEST_MEANANDVARIANCE_H
