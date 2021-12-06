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
#include "PotentialEnergyPivot.h"
#include "BoundSwapPivot.h"
//class ConvexPMF;

class SimplexMCMC: public glp::Simplex {
public:
    typedef BoundSwapPivot Proposal;

//    static constexpr double fractionalK = 0.1;
//    static constexpr double tol = 1e-8;

    class MCMCStatistics {
    public:
        int nSamples = 0;
        int nAccepted = 0;
        int nNonDegenerate = 0;
        int nSwaps = 0;
        int nNulls = 0;
        int nFeasibilityTransitions = 0;

        void update(bool accepted, const Proposal &proposal, bool isSameFeasibilityAsLastLog);
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
    double currentLogProb;
    Proposal proposalFunction;
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


    double lnProb() { return logProbFunc(X()); }

    const std::vector<double> & nextSample();

    using glp::Simplex::pivot;
    void pivot(const Proposal &piv) { pivot(piv.i, piv.j, piv.tableauCol(), piv.leavingVarToUpperBound); }

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

    double infeasibilityGradient(int i) {
        int k = head[i];
        return infeasibilityGradient(beta[i], l[k], u[k]);
    }

    static double infeasibilityGradient(double v, double lowerBound, double upperBound) {
        if(v > upperBound) {
            return 1.0;
        } else if(v < lowerBound) {
            return -1.0;
        }
        return 0.0;
    }


protected:
    bool processProposal(Proposal &proposal);
    Proposal &proposePivot();
    // int proposeColumn();


//    void toCanonicalState();
//    std::vector<int> auxiliaries();

    void updateLPSolution(int j, double deltaj, const glp::FVSVector &pivot);
    void revertLPSolution(int j, double deltaj, const glp::FVSVector &pivot);
    void updateLPSolution(int j, double deltaj, const glp::SparseVec &pivot);
    void revertLPSolution(int j, double deltaj, const glp::SparseVec &pivot);


    double infeasibility();

    void checkLPSolution();

};


#endif //GLPKTEST_SIMPLEXMCMC_H
