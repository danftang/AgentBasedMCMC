//
// Created by daniel on 11/05/2021.
//

#ifndef GLPKTEST_PROPOSALPIVOT_H
#define GLPKTEST_PROPOSALPIVOT_H


#include <utility>
#include <tuple>
#include <vector>
#include <cmath>

#include "glpkpp.h"
#include "constants.h"

class ProposalPivot {
public:

    int i;
    int j;
    bool    leavingVarToUpperBound;
    double  deltaj;                     // the change in the value of this column on pivoting

    glp::Simplex &simplex;
    double logAcceptanceContribution;   // the log probability of this transition minus log-prob of the reverse transition
    std::vector<double> col;            // the j'th column of the tableau
    std::vector<int>    nonZeroRows;    // indices of 'col' that are non-zero


    ProposalPivot(glp::Simplex &simplex): simplex(simplex) { }
    ProposalPivot(glp::Simplex &simplex, int i, int j): i(i), simplex(simplex) {
        if(j>0) setCol(j);
    }

    ProposalPivot(glp::Simplex &simplex, int i, int j, std::vector<double> column):
    i(i), j(j),
    simplex(simplex),
    col(std::move(column)),
    logAcceptanceContribution(0.0),
    nonZeroRows() {
        initNonZeroRows();
    }

    void setCol(int j) {
        this->j = j;
        col = simplex.tableauCol(j);
        initNonZeroRows();
    }


    void initNonZeroRows();

    std::multimap<double,int> getPivotsByDeltaJ();
    std::multimap<double, int> getPivotsByInfeasibility();

    double colInfeasibilityGradient(double deltaj);
    static double infeasibilityGradient(double v, double lowerBound, double upperBound);

    //    operator std::tuple<int &,int &>() { return std::tie(i,j); }
    std::vector<double> reverseCol() const;
//    ProposalPivot reverse() const { return ProposalPivot(i, j, reverseCol()); }
};



#endif //GLPKTEST_PROPOSALPIVOT_H
