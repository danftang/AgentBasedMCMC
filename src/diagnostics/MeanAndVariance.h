// Calculates the mean and varaince of each dimension in a multivaraite object
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


    template<typename DATAPOINTS>
    MeanAndVariance(const DATAPOINTS &datapoints): MeanAndVariance() {
        addDatapoints(datapoints);
    }

    template<typename CONTAINER>
    void addDatapoints(const CONTAINER &datapoints) {
        for(const auto &sample : datapoints) (*this)(sample);
    }

    int dataDimension() const {
        return sum.size();
    }

//    MeanAndVariance(int dimension): nSamples(0), sum(0.0, dimension), sumOfSquares(0.0, dimension) { }

    template<typename T, typename = std::enable_if_t<std::is_convertible_v<decltype(std::declval<T>()[0]),double>>>
    bool operator()(const T &dataPoint) {
        if(nSamples == 0) initArrays(dataPoint.size());
        if(dataPoint.size() != sum.size()) throw("Dimension of dataPoint does not match this logger.");
        for(int i=0; i < dataPoint.size(); ++i) {
            double xi = dataPoint[i];
            sum[i] += xi;
            sumOfSquares[i] += xi*xi;
        }
        ++nSamples;
        return true;
    }

//    template<typename T, typename = std::void_t<decltype(std::declval<MeanAndVariance>()(std::declval<T>()))>>
    auto consumer() {
        return [this](const auto &sample) { return (*this)(sample); };
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

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & nSamples & sum & sumOfSquares;
    }
};


#endif //GLPKTEST_MEANANDVARIANCE_H
