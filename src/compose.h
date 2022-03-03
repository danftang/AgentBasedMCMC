//
// Created by daniel on 30/08/2021.
//

#ifndef GLPKTEST_COMPOSE_H
#define GLPKTEST_COMPOSE_H

#include <functional>

template<typename OUT, typename MID>
std::function<OUT()> operator *(std::function<OUT(const MID &)> g, std::function<MID()> f) {
    return [G = std::move(g), F = std::move(f)]() { return G(F()); };
}

#endif //GLPKTEST_COMPOSE_H
