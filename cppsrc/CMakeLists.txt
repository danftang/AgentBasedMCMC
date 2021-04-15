cmake_minimum_required(VERSION 3.17)
project(glpkTest)

set(CMAKE_CXX_STANDARD 20)
# link_libraries(glpk)


include_directories(glpk glpk/amd glpk/api glpk/bflib  glpk/colamd  glpk/draft  glpk/env  glpk/intopt  glpk/minisat  glpk/misc  glpk/mpl  glpk/npp glpk/proxy glpk/simplex  glpk/zlib)

file(GLOB GLPK_SOURCES
        glpk/*/*.h
        glpk/*/*.c)

add_executable(glpkTest main.cpp ${GLPK_SOURCES} GlpProblem.cpp GlpProblem.h SparseVec.cpp SparseVec.h GlpTableau.cpp GlpTableau.h GlpSimplex.cpp GlpSimplex.h)
