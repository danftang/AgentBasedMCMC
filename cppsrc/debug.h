//
// Created by daniel on 19/07/2021.
//

#ifndef GLPKTEST_DEBUG_H
#define GLPKTEST_DEBUG_H

#ifdef NDEBUG
#define debug(expr) static_cast<void>(0)
#define fail() static_cast<void>(0)
#else
#include "assert.h"

#define debug(expr) expr
#define fail() assert(false)
#endif

#endif //GLPKTEST_DEBUG_H
