//
// Created by daniel on 15/07/2021.
//

#ifndef GLPKTEST_POISSONSTATE_H
#define GLPKTEST_POISSONSTATE_H

#include "ModelState.h"
#include "PMF.h"

template<typename AGENT>
class PoissonState {
public:
    int nSamples;
    ModelState<AGENT> stateCounts;

    PoissonState(): nSamples(0), stateCounts() { }

    PoissonState(std::function<double (const AGENT &)> probability) {
        nSamples = 1e6;
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            AGENT agent(agentId);
            double l = probability(agent);
            if(l > 0.0) stateCounts[agent] = l*nSamples - 0.5;
        }
    }

    PoissonState &operator +=(const std::vector<AGENT> &stateVector) {
        stateCounts += stateVector;
        ++nSamples;
        return *this;
    }


    // Addition for aggregation of states
    PoissonState &operator +=(const ModelState<AGENT> &other) {
        stateCounts += other;
        ++nSamples;
        return *this;
    }


    double lambda(int agentId) const {
        return (stateCounts[agentId]+0.5)/nSamples;
    }


    double logProb(const std::vector<double> &instance) const {
        double logP = 0.0;
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            double k = fabs(instance[agentId]);
            double l = lambda(agentId);
            logP += k*log(l) - l - lgamma(k+1); // log of Poisson
        }
        return logP;
    }

//    std::vector<double> nextSample() { return ((const PoissonState<AGENT> *)this)->nextSample(); }
//    double logProb(const std::vector<double> &X) { return ((const PoissonState<AGENT> *)this)->logProb(X); }

    ModelState<AGENT> nextSample() const {
        ModelState<AGENT> state;
        for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
            state[agentId] = Random::nextPoisson(lambda(agentId));
        }
        return state;
    }

    friend std::ostream &operator <<(std::ostream &out, const PoissonState<AGENT> &poissonState) {
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            if(poissonState.stateCounts[agentId] != 0.0) out << AGENT(agentId) << " -> " << poissonState.lambda(agentId) << " ";
        }
        return  out;
    }


    friend Gnuplot &operator<<(Gnuplot &gp, const PoissonState<AGENT> &state) {
        typedef std::tuple<double, double, double, double, double> HeatRecord;
        std::vector<std::vector<HeatRecord>> heatData;

        for (int x = 0; x < PredPreyAgent::GRIDSIZE; ++x) {
            std::vector<HeatRecord> &record = heatData.emplace_back();
            for (int y = 0; y < PredPreyAgent::GRIDSIZE; ++y) {
                double lPrey = state.lambda(PredPreyAgent(x, y, PredPreyAgent::PREY));
                double lPred = state.lambda(PredPreyAgent(x, y, PredPreyAgent::PREDATOR));
                record.emplace_back(x, y, std::min(lPrey, 1.0) * 200.0, 0.0, std::min(lPred, 1.0) * 200.0);
            }
        }

        gp << "plot [-0.5:" << PredPreyAgent::GRIDSIZE - 0.5 << "][-0.5:" << PredPreyAgent::GRIDSIZE - 0.5 << "] ";
        gp << "'-' with rgbimage notitle\n";
        gp.send2d(heatData);
        return gp;
    }


};

#endif //GLPKTEST_POISSONSTATE_H
