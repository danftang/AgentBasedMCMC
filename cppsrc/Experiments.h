//
// Created by daniel on 06/05/2021.
//

#ifndef GLPKTEST_EXPERIMENTS_H
#define GLPKTEST_EXPERIMENTS_H


class Experiments {
public:
    static void CatMouseExpt();
    static void RandomWalk();

    static double nullPMF(const std::vector<double> &X) { return 0.0; }
};


#endif //GLPKTEST_EXPERIMENTS_H
