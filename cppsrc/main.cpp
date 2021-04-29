#include <iostream>
#include "glpkppinclude/glpkpp.h"

using namespace glpkpp;

int main() {

    GlpProblem myProb;

    myProb.addConstraint(1.0*X(1) + 1.0*X(2) + 1.0*X(3) <= 100.0);
    myProb.addConstraint(10.0*X(1) + 4.0*X(2) + 5.0*X(3) <= 600.0);
    myProb.addConstraint(2.0*X(1) + 2.0*X(2) + 6.0*X(3) <= 300.0);
    myProb.addConstraint(0.0 <= 1.0*X(1));
    myProb.addConstraint(0.0 <= 1.0*X(2));
    myProb.addConstraint(0.0 <= 1.0*X(3));
    myProb.setObjective(10.0*X(1) + 6.0*X(2) + 4.0*X(3));
    myProb.stdBasis();
    myProb.warmUp();

    GlpSimplex mySimplex(myProb);

    std::cout << myProb;
    std::cout << mySimplex << std::endl;

    mySimplex.pivot(3,3,true);

    std::cout << mySimplex << std::endl;

    return 0;
}
