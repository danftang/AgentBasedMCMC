cmake_minimum_required(VERSION 3.17)
project(glpkTest)

set(CMAKE_CXX_STANDARD 17)
link_directories(lib)
link_libraries(glpkpp)

#include_directories(glpk glpk/amd glpk/api glpk/bflib  glpk/colamd  glpk/draft  glpk/env  glpk/intopt  glpk/minisat  glpk/misc  glpk/mpl  glpk/npp glpk/proxy glpk/simplex  glpk/zlib)
include_directories(glpkppinclude)

#file(GLOB GLPK_SOURCES
#        glpk/*/*.h
#        glpk/*/*.c)

#add_executable(glpkTest main.cpp ${GLPK_SOURCES} GlpProblem.cpp GlpProblem.h SparseVec.cpp SparseVec.h GlpTableau.cpp GlpTableau.h GlpSimplex.cpp GlpSimplex.h SimplexMCMC.cpp SimplexMCMC.h)
add_executable(glpkTest main.cpp Agent.h ABM.h agents/CatMouseAgent.cpp agents/CatMouseAgent.h Observation.h ABMProblem.h Event.h Trajectory.h State.h Experiments.cpp Experiments.h SimplexMCMC.cpp StlStream.h BasisProbability.cpp BasisProbability.h ProposalPivot.h Random.cpp Random.h ProposalPivot.cpp StateTrajectory.h ColumnPivot.h ColumnPivot.cpp agents/PredPreyAgent.cpp agents/PredPreyAgent.h ModelState.h ProbabilisticColumnPivot.cpp ProbabilisticColumnPivot.h)
