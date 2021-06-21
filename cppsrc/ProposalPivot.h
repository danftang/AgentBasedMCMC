//
// Created by daniel on 11/05/2021.
//

#ifndef GLPKTEST_PROPOSALPIVOT_H
#define GLPKTEST_PROPOSALPIVOT_H


#include <utility>
#include <tuple>
#include <vector>

#include "glpkpp.h"

class ProposalPivot {
public:

    int i;
    int j;
    bool leavingVarToUpperBound;
    std::vector<double> col;        // the j'th column of the tableau
    double logAcceptanceContribution;  // the log probability of proposing this transition over probability of proposing reverse transition
    double              deltaj;      // the change in the value of this column on pivoting

    ProposalPivot() = default;
    ProposalPivot(int i, int j, std::vector<double> column): i(i), j(j), col(std::move(column)), logAcceptanceContribution(0.0) { }

    //    operator std::tuple<int &,int &>() { return std::tie(i,j); }
    std::vector<double> reverseCol() const;
//    ProposalPivot reverse() const { return ProposalPivot(i, j, reverseCol()); }
};



#endif //GLPKTEST_PROPOSALPIVOT_H
