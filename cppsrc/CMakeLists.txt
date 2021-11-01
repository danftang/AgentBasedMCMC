cmake_minimum_required(VERSION 3.17)
project(glpkTest)

set(CMAKE_CXX_STANDARD 17)
link_directories(lib)
link_libraries(glpkpp boost_system boost_iostreams pthread)

#include_directories(glpk glpk/amd glpk/api glpk/bflib  glpk/colamd  glpk/draft  glpk/env  glpk/intopt  glpk/minisat  glpk/misc  glpk/mpl  glpk/npp glpk/proxy glpk/simplex  glpk/zlib)
include_directories(glpkppinclude)

#file(GLOB GLPK_SOURCES
#        glpk/*/*.h
#        glpk/*/*.c)

#add_executable(glpkTest main.cpp ${GLPK_SOURCES} GlpProblem.cpp GlpProblem.h SparseVec.cpp SparseVec.h GlpTableau.cpp GlpTableau.h GlpSimplex.cpp GlpSimplex.h SimplexMCMC.cpp SimplexMCMC.h)
add_executable(glpkTest main.cpp Agent.h ABM.h agents/CatMouseAgent.cpp agents/CatMouseAgent.h Observation.h ABMProblem.h Event.h Trajectory.h State.h Experiments.cpp Experiments.h SimplexMCMC.cpp StlStream.h ProposalPivot.h Random.cpp Random.h ProposalPivot.cpp StateTrajectory.h Phase2Pivot.h Phase2Pivot.cpp agents/PredPreyAgent.cpp agents/PredPreyAgent.h ModelState.h constants.h Phase1Pivot.cpp Phase1Pivot.h PotentialEnergyPivot.cpp PotentialEnergyPivot.h DataAssimilation.h PoissonState.h debug.h AssimilationWindow.h BinomialObservation.h ActFermionicDistribution.h UnitTests.h ConvexPMF.h TrajectoryPriorDistribution.h agents/BinomialAgent.cpp agents/BinomialAgent.h DeltaPMF.h ConvexPolyhedron.h ABMConstraints.h PoissonPMF.h BinomialDistribution.h AssimilationProblem.h ABMImportanceSampler.h RejectionSampler.h ConvexPMFProduct.cpp ConvexPMFProduct.h AgentStateObservation.h Distribution.h PMFProduct.h BinarySolutionSet.h ExpectationProblem.h SampleStatistics.h TrajectoryLikelihoodPMF.h TrajectorySampler.h ExactSolver.h MCMCSolver.h Plotter.h IntSampleStatistics.h ABMImportanceSampler.cpp ABMEndStateSolver.h ConvexPMFBase.h BernoulliModelState.h MCMCSampler.h compose.h ModelStateSampleStatistics.h Sampler.h diagnostics/KLDivergence.h diagnostics/Autocorrelation.h diagnostics/MeanAndVariance.h diagnostics/ValarrayLogger.h diagnostics/ScaleReduction.h GnuplotExtensions.h diagnostics/Dataflow.h)
