//
// Created by daniel on 23/07/2021.
//

#ifndef GLPKTEST_ACTFERMIONICDISTRIBUTION_H
#define GLPKTEST_ACTFERMIONICDISTRIBUTION_H


#include <vector>
#include <cassert>
#include "Random.h"
#include "StlStream.h"

// This is the distribution of acts of a number of act-fermionic agents
// in the same state. Given m agents and a vector of act probabilities pi=<p_1...p_N>
// then in the non-Fermionic case, the probability of an m-dimensional vector of acts
// for each of the m agents, in order, A=<a_1...a_m> is given by
//
// P(A|pi) = \prod_{i=1}^m p_{a_i}
//
// In the Fermionic case, we assume that all act occupations are at most 1, so we
// restrict ourselves to the set of vectors whose elements are all different,
// which we call C^N_m. So for Fermionic agents we have
//
//                        { (1/S^N_m) \prod_{i=1}^m p_{a_i}     if B \in C^N_m
// P^N_m(B|B \in C^N_m) = {
//                        { 0                               otherwise
// where
// S^N_m = \sum_{X \in C^N_m} \prod_{i=1}^m p_{x_i}
//
// In order to sample from this distribution, given the total number of acts
// that must be 'active', we note that the marginalised probability that the last
// act is inactive is given by
//
// P(b_N=0) = S^(N-1)_m/S^N_m
//
// So we can choose the status of the N'th act with that probbility then if we choose b_N=0 then
// choose the other bits by solving the sub-problem N->N-1 and m->m, or if we choose b_N=1 then
// choose the other bits by solving the sub-problem N->N-1, m->m-1.
//
class ActFermionicDistribution {
public:
    std::vector<std::vector<double>>    s;      // s[i][j] is sum of terms over the first N=i+j+1 acts with m=j-1 agents
    std::vector<double>                 p;      // probabilities of each of N outcomes

    ActFermionicDistribution(std::vector<double> probabilities): p(probabilities) {
        s.push_back(std::vector<double>({p[0]}));
//        std::cout << "Distribution = " << p << std::endl;
    }


    int nActs() { return p.size(); }

    // log probability
    double operator()(const std::vector<bool> &X, int m=-1) {
        if(m == -1) m = std::count_if(X.begin(), X.end(), [](bool b) { return b; });
        double P = 1.0;
        for(int i=0; i<X.size(); ++i) {
            if(X[i]) P *= p[i];
        }
        return P/sigma(X.size(), m);
    }

    // return a bit-vector <b_1...b_N> with m bits true sampled from this distribution
    // or return the empty vector if no such vector exists (possible if some p[i] are 0.0)
    std::vector<bool> sampleUnordered(int m) {
        if(m > p.size()) return std::vector<bool>();
        std::vector<bool> actOccupation(p.size(), false);
        int N=p.size();
        while(m>0) {
//            double pUnoccupied = sigma(act, m)/sigma(act+1, m);
            double sNm = sigma(N,m);
            if(sNm == 0.0) return std::vector<bool>();
            if(N == m || Random::nextDouble() > sigma(N-1, m)/sNm) {
                actOccupation[N-1] = true;
                --m;
            }
            --N;
        }
        return actOccupation;
    }


    double sigma(int N, int m) {
        assert(N >= 1 && m >=1);
        assert(N >= m);
        int i = N-m;
        int j = m-1;
        if(i >= s.size() || j >= s[i].size()) expandS(i,j);
//        std::cout << "Sigma " << N << ", " << m << " = " << s[i][j] << std::endl;
        return s[i][j];
    }



    void expandS(int iMax, int jMax) {
        if(iMax >= s.size()) s.resize(iMax+1);
        for(int i=0; i<=iMax; ++i) {
            std::vector<double> &si = s[i];
            if(si.size() <= jMax) {
                if(si.size() == 0) si.push_back(p[i] + s[i-1][0]);
                for(int j=si.size(); j <= jMax; ++j) {
                    double sij = (j+1)*p[i+j]*s[i][j-1];
                    if(i>0) sij += s[i-1][j];
                    si.push_back(sij);
                }
            }
        }
    }
};


#endif //GLPKTEST_ACTFERMIONICDISTRIBUTION_H
