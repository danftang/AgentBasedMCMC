//
// Created by daniel on 18/10/2021.
//

#ifndef GLPKTEST_KLDIVERGENCE_H
#define GLPKTEST_KLDIVERGENCE_H

#include <vector>
#include <deque>
#include <valarray>
#include "../gnuplot-iostream/gnuplot-iostream.h"

/////////////////////////////////////////////////////////////
// Calculates the KL divergence D(Samples|Real distribution)
// between the normalised sum of deltas on the sample points and the
// target distribution. This is calculable (in this direction)
// since it is just
//
// D(s_0....s_n|P) = sum_{i=0}^n -log(P(s_i))/n
//
// Since we don't want to store all samples, and we're only
// interested in the sample point probabilities, just pipe the sample
// probabilities in to an instance of this class. The class will
// calculate divergence in a sliding window, logging the value
// at every m samples. The window length and log
// period are supplied to the constructor.
/////////////////////////////////////////////////////////////
class KLDivergence {
public:
    const int               logPeriod;
    std::vector<double>     divergenceLog;
    std::valarray<double>   window;
    int                     windowStart;
    int                     timeToNextLog;

    KLDivergence(int windowLength, int logPeriod):
        logPeriod(logPeriod),
        window(0.0, windowLength),
        windowStart(0),
        timeToNextLog(windowLength) {
        std::cout << "Creating KL-Divergence logger with window length " << windowLength << "and period " << logPeriod << std::endl;
    }


    KLDivergence &operator <<(double sampleLogProb) {
        window[windowStart] = sampleLogProb;
        windowStart = (windowStart + 1)%window.size();
        if(--timeToNextLog == 0) {
            divergenceLog.push_back(-window.sum()/window.size());
            timeToNextLog = logPeriod;
        }
        return *this;
    }


    void plot() {
        Gnuplot gp;
        gp << "plot '-' with lines\n";
        gp.send1d(divergenceLog);
    }
};


#endif //GLPKTEST_KLDIVERGENCE_H
