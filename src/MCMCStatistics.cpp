//
// Created by daniel on 30/03/2022.
//

#include <iostream>
#include "MCMCStatistics.h"

MCMCStatistics::MCMCStatistics() {
    reset();
}

void MCMCStatistics::addSample(bool proposalAccepted, bool startStateIsFeasible, bool proposalIsFeasible) {
    ++nProposals[startStateIsFeasible][proposalIsFeasible];
    if(proposalAccepted) {
        ++nAccepted[startStateIsFeasible][proposalIsFeasible];
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

int MCMCStatistics::totalProposals() const {
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
    out << "Total feasible samples      " << stats.nFeasibleSamples() << std::endl;
    out << "Proportion infeasible       " << stats.nInfeasibleSamples()*100.0/ stats.totalProposals() << "%" << std::endl;
    out << "Mean feasible run-length    " << stats.nFeasibleSamples()/(stats.nAccepted[false][true]+1.0) << std::endl;
    out << "Mean infeasible run-length  " << stats.nInfeasibleSamples()*1.0/stats.nAccepted[true][false] << std::endl;
//    out << "Mean acceptance             " << stats.totalAccepted()*100.0/ stats.totalProposals() << "%" << std::endl;
//    out << "Proposals" << std::endl;
//    out << "   total                    " << stats.totalProposals() << std::endl;
//    out << "   feasible-feasible        " << stats.nProposals[true][true]*100.0/ stats.totalProposals() << "%" << std::endl;
//    out << "   feasible-infeasible      " << stats.nProposals[true][false]*100.0/ stats.totalProposals() << "%" << std::endl;
//    out << "   infeasible-infeasible    " << stats.nProposals[false][false]*100.0/ stats.totalProposals() << "%" << std::endl;
//    out << "   infeasible-feasible      " << stats.nProposals[false][true]*100.0/ stats.totalProposals() << "%" << std::endl;
    out << "Acceptance" << std::endl;
    out << "   mean                     " << stats.totalAccepted()*100.0/ stats.totalProposals() << "%" << std::endl;
    out << "   feasible-feasible        " << stats.nAccepted[true][true]*100.0/stats.nProposals[true][true] << "%" << std::endl;
    out << "   feasible-infeasible      " << stats.nAccepted[true][false]*100.0/stats.nProposals[true][false] << "%" << std::endl;
    out << "   infeasible-infeasible    " << stats.nAccepted[false][false]*100.0/stats.nProposals[false][false] << "%" << std::endl;
    out << "   infeasible-feasible      " << stats.nAccepted[false][true]*100.0/stats.nProposals[false][true] << "%" << std::endl;

    return out;
}
