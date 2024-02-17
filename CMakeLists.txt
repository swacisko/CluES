cmake_minimum_required(VERSION 3.14)
project(CluES_exp)

set(CMAKE_CXX_STANDARD 20)

#******************************************************************************************  Google-tests
#configure_file(CMakeLists.txt.gtest googletest-download/CMakeLists.txt)
#execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#        RESULT_VARIABLE result
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#if(result)
#    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
#endif()
#execute_process(COMMAND ${CMAKE_COMMAND} --build .
#        RESULT_VARIABLE result
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#if(result)
#    message(FATAL_ERROR "Build step for googletest failed: ${result}")
#endif()

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)
#******************************************************************************************  Google-tests


##********************************************************************* google-benchmark
#set(CMAKE_BUILD_TYPE Release) # this enables building googlebenchmark in release mode

#configure_file(CMakeLists.txt.gbench googlebenchmark-download/CMakeLists.txt)
#execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . # build in debug mode
##execute_process(COMMAND ${CMAKE_COMMAND} "${CMAKE_GENERATOR}" . # build in release mode
#        RESULT_VARIABLE result
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googlebenchmark-download )
#if(result)
#    message(FATAL_ERROR "CMake step for googlebenchmark failed: ${result}")
#endif()
#execute_process(COMMAND ${CMAKE_COMMAND} --build .
#        RESULT_VARIABLE result
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googlebenchmark-download )
#if(result)
#    message(FATAL_ERROR "Build step for googlebenchmark failed: ${result}")
#endif()
#
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googlebenchmark-src
        ${CMAKE_CURRENT_BINARY_DIR}/googlebenchmark-build
        EXCLUDE_FROM_ALL)


#********************************************************************* google-benchmark


#SET(CMAKE_CXX_FLAGS "-std=c++20 -O0")
SET(CMAKE_CXX_FLAGS "-std=c++20 -O3 -march=native")

#***************************** Cluster Editing

include_directories(
        include
        include/graphs
        include/combinatorics
        include/utils
        include/datastructures
        include/clues
        include/clues/*
        include/clues/*/*
        include/clues/*/*/*
        include/combinatorics/*
        include/datastructures/*
        experiments
        )

file(GLOB SOURCES
        src/graphs/GraphReader.cpp
        src/graphs/GraphWriter.cpp
        src/graphs/generators/GraphGenerator.cpp
        src/graphs/GraphUtils.cpp
        src/graphs/GraphInducer.cpp
        src/combinatorics/CombinatoricUtils.cpp

        src/datastructures/FAU.cpp

        src/utils/RandomNumberGenerators.cpp
        src/utils/TimeMeasurer.cpp

        src/*
        src/*/*
        src/*/*/*
        src/*/*/*/*

        src/clues/CETestGraphGenerator.cpp
        experiments/bm_exp.cpp
        experiments/Stopwatch.cpp

        ) # include all source files

#***************************** Cluster Editing



add_executable(CluES_exp ${SOURCES})
target_link_libraries(CluES_exp benchmark::benchmark)
target_link_libraries(CluES_exp gtest)
target_link_libraries(CluES_exp gmock)
target_link_libraries(CluES_exp stdc++fs)
