//
// Created by daniel on 20/10/2021.
//

#ifndef GLPKTEST_MEANANDVARIANCE_H
#define GLPKTEST_MEANANDVARIANCE_H


#include <valarray>

class MeanAndVariance {
public:
    int nSamples;
    std::valarray<double> sum;
    std::valarray<double> sumOfSquares;

    MeanAndVariance(int dimension): nSamples(0), sum(0.0, dimension), sumOfSquares(0.0, dimension) { }

    // ARRAY can be any object with .size() and operator[]
    template<typename ARRAY>
    void operator()(const ARRAY &sample) {
        if(sample.size() != sum.size()) throw("Dimension of sample does not match this logger.");
        for(int i=0; i<sample.size(); ++i) {
            double xi = sample[i];
            sum[i] += xi;
            sumOfSquares[i] += xi*xi;
        }
        ++nSamples;
    }

    std::valarray<double> mean() const {
        return sum * (1.0/nSamples);
    }


    std::valarray<double> sampleVariance() const {
        return (sumOfSquares - sum*sum * (1.0/nSamples))*(1.0/(nSamples-1));
    }
};


#endif //GLPKTEST_MEANANDVARIANCE_H
