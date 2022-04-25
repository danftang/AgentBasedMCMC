// Represents various statistics on the performance of a Markov Chain
//
// Created by daniel on 30/03/2022.
//

#ifndef ABMCMC_MCMCSTATISTICS_H
#define ABMCMC_MCMCSTATISTICS_H

#include <ostream>
#include "boost/serialization/access.hpp"

class MCMCStatistics {
public:
    int nAccepted[2][2];    // number of proposals accepted by start feasibility and end feasibility
    int nProposals[2][2];   // total number of proposals by start feasibility and end feasibility
//    double sumOfLogImportances;

    MCMCStatistics();
    void addSample(bool proposalAccepted, bool isCurrentlyFeasible, bool proposalIsFeasible);
    void reset();

    int nSamples() const;
    int totalAccepted() const;
    int nRejected(bool startIsFeasible, bool proposalIsFeasible) const;
    int nInfeasibleSamples() const;
    int nFeasibleSamples() const;

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & nAccepted & nProposals;
    }

};

std::ostream &operator <<(std::ostream &out, const MCMCStatistics &stats);

#endif //ABMCMC_MCMCSTATISTICS_H
