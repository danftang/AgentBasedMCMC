//
// Created by daniel on 09/11/2021.
//

#include <string>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "FiguresForPaper.h"
#include "PredPreyProblem.h"


void FiguresForPaper::generateProblemFile() {
    constexpr int GRIDSIZE = 8;
    constexpr int nTimesteps = 2;
    constexpr double pPredator = 0.03;
    constexpr double pPrey = 0.06;
    constexpr double pMakeObservation = 0.04;
    constexpr double pObserveIfPresent = 0.9;
    std::string filename = "PredPrey" + std::to_string(GRIDSIZE) + "." + std::to_string(nTimesteps) + ".prob";
    std::ofstream file(filename);
    boost::archive::binary_oarchive boostArchive(file);

    PredPreyProblem<GRIDSIZE> problem(nTimesteps, pPredator, pPrey, pMakeObservation, pObserveIfPresent);
    std::cout << problem << std::endl;
    boostArchive << problem;
}

