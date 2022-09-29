cmake_minimum_required(VERSION 3.10)
project(ABMCMC)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
link_directories(lib)
link_libraries(boost_system boost_iostreams boost_filesystem boost_serialization pthread)

add_executable(abmcmc main.cpp MCMCStatistics.cpp)
#add_executable(abmcmc main.cpp MCMCStatistics.cpp Experiments.cpp agents/CatMouseAgent.cpp)
