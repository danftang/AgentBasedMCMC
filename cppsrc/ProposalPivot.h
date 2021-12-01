//
// Created by daniel on 11/05/2021.
//

#ifndef GLPKTEST_PROPOSALPIVOT_H
#define GLPKTEST_PROPOSALPIVOT_H


#include <utility>
#include <tuple>
#include <vector>
#include <cmath>
#include <map>
#include <iostream>

#include "constants.h"
class SimplexMCMC;

class ProposalPivot {
public:

    int i;
    int j;
    bool    leavingVarToUpperBound;
    double  deltaj;                     // the change in the value of this column on pivoting

    SimplexMCMC &simplex;
    double logAcceptanceContribution;   // the log probability of this transition minus log-prob of the reverse transition
    std::vector<double> col;            // the j'th column of the tableau
    std::vector<int>    nonZeroRows;    // indices of 'col' that are non-zero


    ProposalPivot(SimplexMCMC &simplex): simplex(simplex) { }
    ProposalPivot(SimplexMCMC &simplex, int i, int j): i(i), simplex(simplex) {
        if(j>0) {
            setCol(j);
        } else {
            this->j = j;
        }
    }

    ProposalPivot(SimplexMCMC &simplex, int i, int j, std::vector<double> column):
    i(i), j(j),
    simplex(simplex),
    col(std::move(column)),
    logAcceptanceContribution(0.0),
    nonZeroRows() {
        initNonZeroRows();
    }

    void setCol(int j);
    void initNonZeroRows();
    void clearCol();

    std::multimap<double,int> getPivotsByDeltaJ();
    std::multimap<double, int> getPivotsByInfeasibility();

    double colInfeasibilityGradient(double deltaj);
    static double infeasibilityGradient(double v, double lowerBound, double upperBound) {
        if(v > upperBound) {
            return 1.0;
        } else if(v < lowerBound) {
            return -1.0;
        }
        return 0.0;
    }

    //    operator std::tuple<int &,int &>() { return std::tie(i,j); }
    std::vector<double> reverseCol() const;
//    ProposalPivot reverse() const { return ProposalPivot(i, j, reverseCol()); }

protected:

    bool isActive(int pmfIndex);
    void setToPivotIndex(int pivotIndex);

    std::vector<double> infeasibilityCost();
};



#endif //GLPKTEST_PROPOSALPIVOT_H
