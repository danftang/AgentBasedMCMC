//
// Created by daniel on 21/10/2021.
//

#ifndef GLPKTEST_SCALEREDUCTION_H
#define GLPKTEST_SCALEREDUCTION_H

/////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the potential scale reduction as defined in
// Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
// Chapter 11, page 285.
//
// R^hat =  B/(nW) + (n-1)/n
//
// where n is the number of samples per chian, B is the between-chains variance (i.e. n times the
// vaiance of the means of each chain) and W is the within chains variance (i.e. the mean of the
// sample variances of each chain).
//
// Input is a vector of chain MeanAndVariances, output is a vector whose elements are the
// potential scale reductions in each dimension of the samples.
/////////////////////////////////////////////////////////////////////////////////////////////
std::valarray<double> gelmanScaleReduction(const std::vector<MeanAndVariance> &sampleStats) {
    if(sampleStats.size() == 0) return std::valarray<double>();
    int n = sampleStats.front().nSamples;
    int m = sampleStats.size(); // number of chains
    int d = sampleStats.front().sum.size(); // dimension of a sample;
    std::valarray<double> allChainsMean(0.0, d);
    std::valarray<double> W(0.0, d); // within chains variance;
    std::valarray<double> B(0.0, d); // between chains variance;
    //std::valarray<double> allChainsMean(0.0, d);

    for(auto &chain: sampleStats) {
        allChainsMean += chain.mean();
        W += chain.sampleVariance();
    }
    allChainsMean *= 1.0/m;
    W *= 1.0/m;

    for(auto &chain: sampleStats) B += pow(chain.mean() - allChainsMean,2);
    B *= n/(m-1.0);

    return B / (W*(1.0*n) + 1e-15) + ((n-1.0)/n);
}

#endif //GLPKTEST_SCALEREDUCTION_H
