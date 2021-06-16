#include <iostream>
#include "glpkppinclude/glpkpp.h"
#include "agents/CatMouseAgent.h"
#include "ABMProblem.h"
#include "Trajectory.h"
#include "Experiments.h"
#include "StlStream.h"
#include "Random.h"

using glp::X;



int main() {
//    std::multimap<double, int> myMap();
//
//
//
//    myMap.emplace(0.01,1);
//    myMap.emplace(-0.01,2);
//    myMap.emplace(0.01,3);
//    myMap.emplace(-0.01,4);
//
//    for(auto [key,val] : myMap) {
//        std::cout << key << "->" << val << std::endl;
//    }
//
//    auto lb = myMap.lower_bound(0.0);
//    std::cout << lb->second << std::endl;
//    std::cout << -0.0 << " " << 0.0 << std::endl;

    Experiments::PredPreyExpt();
//    Experiments::CatMouseExpt();
//    Experiments::RandomWalk();

    return 0;
}

