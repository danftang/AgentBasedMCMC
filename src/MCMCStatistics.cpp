//
// Created by daniel on 30/03/2022.
//

#include <iostream>
#include "MCMCStatistics.h"

MCMCStatistics::MCMCStatistics() {
    reset();
}

void MCMCStatistics::addSample(bool proposalAccepted, bool isCurrentlyFeasible, bool proposalIsFeasible) {
    ++nProposals[isCurrentlyFeasible][proposalIsFeasible];
    if(proposalAccepted) {
        ++nAccepted[isCurrentlyFeasible][proposalIsFeasible];
    }
}

void MCMCStatistics::reset() {
    for(int start = 0; start<2; ++start) {
        for(int end = 0; end<2; ++end) {
            nProposals[start][end] = 0;
            nAccepted[start][end] = 0;
        }
    }
}

int MCMCStatistics::nSamples() const {
    return nProposals[true][true] + nProposals[true][false] + nProposals[false][true] + nProposals[false][false];
}

int MCMCStatistics::nRejected(bool startIsFeasible, bool proposalIsFeasible) const {
    return nProposals[startIsFeasible][proposalIsFeasible] - nAccepted[startIsFeasible][proposalIsFeasible];
}

int MCMCStatistics::totalAccepted() const {
    return nAccepted[true][true] + nAccepted[true][false] + nAccepted[false][true] + nAccepted[false][false];
}

int MCMCStatistics::nInfeasibleSamples() const {
    return nProposals[false][false] + nAccepted[true][false] + nRejected(false,true);
}

int MCMCStatistics::nFeasibleSamples() const {
    return nProposals[true][true] + nAccepted[false][true] + nRejected(true,false);
}

std::ostream &operator <<(std::ostream &out, const MCMCStatistics &stats) {
    out << "Total samples                    " << stats.nSamples() << std::endl;
    out << "Total feasible samples           " << stats.nFeasibleSamples() << std::endl;
//    out << "Exp-log-mean importance          " << stats.meanImportance() <<  std::endl;
    out << "proportion infeasible            " << stats.nInfeasibleSamples()*100.0/stats.nSamples() << "%" << std::endl;
    out << "mean feasible run-length         " << stats.nFeasibleSamples()/(stats.nAccepted[false][true]+1.0) << std::endl;
    out << "mean infeasible run-length       " << stats.nInfeasibleSamples()*1.0/stats.nAccepted[true][false] << std::endl;
    out << "mean acceptance                  " << stats.totalAccepted()*100.0/stats.nSamples() << "%" << std::endl;
    out << "   feasible-feasible             " << stats.nAccepted[true][true]*100.0/stats.nProposals[true][true] << "%" << std::endl;
    out << "   feasible-infeasible           " << stats.nAccepted[true][false]*100.0/stats.nProposals[true][false] << "%" << std::endl;
    out << "   infeasible-infeasible         " << stats.nAccepted[false][false]*100.0/stats.nProposals[false][false] << "%" << std::endl;
    out << "   infeasible-feasible           " << stats.nAccepted[false][true]*100.0/stats.nProposals[false][true] << "%" << std::endl;
    return out;
}
