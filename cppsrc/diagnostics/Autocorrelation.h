//
// Created by daniel on 20/10/2021.
//

#ifndef GLPKTEST_AUTOCORRELATION_H
#define GLPKTEST_AUTOCORRELATION_H

#include <vector>
#include <valarray>

////////////////////////////////////////////////////////////////////////////
// Calculates the autocovariance, g_t, of discrete samples x_0...x_n as defined in
// Geyer, C.J., 1992. Practical markov chain monte carlo. Statistical science, pp.473-483.
//
// g_t = 1/(n+1) sum_{i=0}^{n-t} (x_i - mu)(x_{i+t} - mu)
// where
// mu = 1/(n+1) sum_{i=0}^n x_i
//
// The lag, t, is calculated at "nLags" different values
// evently spaced from 0 until "maxLagProportion"*(n+1).
//
// Since we're only interested in estimating the overall rate of decay,
// we don't need to calculate g_t at all values of t, so although this
// could be done by using FFT on the zero-padded samples to convert into
// and back from the frequency domain, the saving is negligible for the sizes
// we're interested in (nLags approx 100) so we just calculate explicitly.
//
// The SAMPLE type can be any type that has arithmetic operations, notably
// another valarray.
//
////////////////////////////////////////////////////////////////////////////
template<typename SAMPLE>
std::valarray<SAMPLE> geyerAutocorrelation(const std::valarray<SAMPLE> &samples, int nLags=100, double maxLagProportion=0.9) {
    std::valarray<SAMPLE> autocov(0.0 * samples[0], nLags);
    SAMPLE meanSample = samples.sum()*(1.0/samples.size());
    std::valarray<SAMPLE> Dx = samples - meanSample;

    int stride = (maxLagProportion*samples.size())/nLags;
    for(int j = 0; j < autocov.size(); ++j) {
        int t = j*stride;
        for(int i=0; i+t < samples.size(); ++i) {
            autocov[j] += Dx[i]*Dx[i+t];
        }
        autocov[j] /= samples.size();
    }
    return autocov;
}




#endif //GLPKTEST_AUTOCORRELATION_H