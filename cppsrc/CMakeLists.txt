cmake_minimum_required(VERSION 3.17)
project(glpkTest)

set(CMAKE_CXX_STANDARD 17)
link_directories(lib)
link_libraries(glpkpp boost_system boost_iostreams)

#include_directories(glpk glpk/amd glpk/api glpk/bflib  glpk/colamd  glpk/draft  glpk/env  glpk/intopt  glpk/minisat  glpk/misc  glpk/mpl  glpk/npp glpk/proxy glpk/simplex  glpk/zlib)
include_directories(glpkppinclude)

#file(GLOB GLPK_SOURCES
#        glpk/*/*.h
#        glpk/*/*.c)

#add_executable(glpkTest main.cpp ${GLPK_SOURCES} GlpProblem.cpp GlpProblem.h SparseVec.cpp SparseVec.h GlpTableau.cpp GlpTableau.h GlpSimplex.cpp GlpSimplex.h SimplexMCMC.cpp SimplexMCMC.h)
add_executable(glpkTest main.cpp Agent.h ABM.h agents/CatMouseAgent.cpp agents/CatMouseAgent.h Observation.h ABMProblem.h Event.h Trajectory.h State.h Experiments.cpp Experiments.h SimplexMCMC.cpp StlStream.h ProposalPivot.h Random.cpp Random.h ProposalPivot.cpp StateTrajectory.h Phase2Pivot.h Phase2Pivot.cpp agents/PredPreyAgent.cpp agents/PredPreyAgent.h ModelState.h ProbabilisticColumnPivot.cpp ProbabilisticColumnPivot.h constants.h Phase1Pivot.cpp Phase1Pivot.h PotentialEnergyPivot.cpp PotentialEnergyPivot.h DataAssimilation.h PoissonState.h debug.h AssimilationWindow.h BinomialObservation.h ActFermionicDistribution.h UnitTests.cpp UnitTests.h ConvexPMF.cpp ConvexPMF.h PMF.h)
