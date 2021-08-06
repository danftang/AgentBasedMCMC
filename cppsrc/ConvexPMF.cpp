//
// Created by daniel on 26/07/2021.
//

#include "ConvexPMF.h"

const ConvexPMF ConvexPMF::INVALID_PMF = ConvexPMF([](const std::vector<double> &X) { return NAN; }, 0);
