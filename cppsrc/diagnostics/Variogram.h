//
// Created by daniel on 02/11/2021.
//

#ifndef GLPKTEST_VARIOGRAM_H
#define GLPKTEST_VARIOGRAM_H

////////////////////////////////////////////////////////////////////////////
// Calculates the variogram, V_t, of discrete samples x_0...x_n as defined in
// Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
// Chapter 11, page 286.
//
// V_t = 1/(n-t) sum_{i=0}^{n-t-1} (x_i - x_{i+t})^2
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
std::valarray<SAMPLE> variogram(const std::valarray<SAMPLE> &samples, int nLags=100, double maxLagProportion=0.5) {
    std::valarray<SAMPLE> vario(0.0 * samples[0], nLags);

    int stride = (maxLagProportion*samples.size())/(nLags-1.0);
    for(int j = 0; j < nLags; ++j) {
        int t = j*stride;
        for(int i=0; i < samples.size()-t; ++i) {
            auto dxt = samples[i] - samples[i+t];
            vario[j] += dxt*dxt;
        }
        vario[j] /= samples.size()-t;
    }
    return vario;
}


#endif //GLPKTEST_VARIOGRAM_H
