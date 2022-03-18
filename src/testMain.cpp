//
// Created by daniel on 18/03/2022.
//

#include <map>
#include <iostream>
#include "Constraint.h"

int main(int argc, char *argv[]) {
    Constraint<int> constraint = (3 <= 4*X(1) + 5*X(2) <= 6);
    std::cout << constraint << std::endl;

}