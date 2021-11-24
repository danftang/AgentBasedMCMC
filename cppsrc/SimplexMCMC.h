//
// Created by daniel on 15/04/2021.
//

#ifndef GLPKTEST_SIMPLEXMCMC_H
#define GLPKTEST_SIMPLEXMCMC_H

#include <vector>
#include <set>
#include <map>
#include <random>
#include <cassert>
#include "glpkpp.h"
#include "ProposalPivot.h"
#include "ConvexPMF.h"
#include "NormalisedPotentialEnergyPivot.h"
//class ConvexPMF;

class SimplexMCMC: public glp::Simplex {
public:
//    static constexpr double fractionalK = 0.1;
//    static constexpr double tol = 1e-8;

    using glp::Simplex::nNonBasic;

    class MCMCStatistics {
    public:
        int nSamples = 0;
        int nAccepted = 0;
        int nNonDegenerate = 0;
        int nSwaps = 0;
        int nNulls = 0;
        int nFeasibilityTransitions = 0;

        void update(bool accepted, const ProposalPivot &proposal, bool isSameFeasibilityAsLastLog);
        friend std::ostream &operator <<(std::ostream &out, const MCMCStatistics &stats);
    private:
        friend class boost::serialization::access;

        template <typename Archive>
        void serialize(Archive &ar, const unsigned int version) {
            ar & nSamples & nAccepted & nNonDegenerate & nSwaps & nNulls & nFeasibilityTransitions;
        }

    };

    MCMCStatistics feasibleStatistics;
    MCMCStatistics infeasibleStatistics;

    std::function<double (const std::vector<double> &)> logProbFunc;
    NormalisedPotentialEnergyPivot proposalFunction;
    bool lastSampleWasAccepted = true;
//    BasisProbability probability;


    SimplexMCMC(const glp::Problem &prob,
                std::function<double(const std::vector<double> &)> logProb,
                const std::vector<double> &initialState = std::vector<double>());

//    template<typename DOMAIN>
//    SimplexMCMC(const ConvexPMF<DOMAIN> &pmf, const DOMAIN &initialState)
//    : SimplexMCMC(
//            pmf.convexSupport.toLPProblem(),
//            [logP = pmf.extendedLogProb](const std::vector<double> &X) { return logP(reinterpret_cast<const DOMAIN &>(X)); },
//            initialState
//    ) { }


    double lnDegeneracyProb();
    double lnProb() { return logProbFunc(X()); }
    // double lnFractionalPenalty();

    const std::vector<double> & nextSample();
//    double reverseTransitionProb(ProposalPivot proposal);
    using glp::Simplex::pivot;
    void pivot(const ProposalPivot &piv) {
        this->glp::Simplex::pivot(piv.i, piv.j, piv.col, piv.leavingVarToUpperBound);
    }


//    std::vector<int> calcPivotRows(int j, const std::vector<double> &colVec);
//    std::vector<int> calcPivotRows(int j) { return calcPivotRows(j,tableauCol(j)); }

    void randomWalk();

    static glp::Problem &initialiseProblem(glp::Problem &lp);

    void findFeasibleStartPoint(); // phase 1

    void setLPState(const std::vector<double> &lpState);

    // TEST STUFF
    int countFractionalPivCols();
    int infeasibilityCount();
    bool abmSanityChecks();

    bool solutionIsPrimaryFeasible();
    bool solutionIsInteger();
protected:
    bool processProposal(const ProposalPivot &proposal);
    const ProposalPivot &proposePivot();
    // int proposeColumn();


//    void toCanonicalState();
//    std::vector<int> auxiliaries();

    void updateLPSolution(const ProposalPivot &pivot);
    void revertLPSolution(const ProposalPivot &pivot);


    double infeasibility();

};


#endif //GLPKTEST_SIMPLEXMCMC_H
