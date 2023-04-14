cmake_minimum_required(VERSION 3.14)
project(CluES)

set(CMAKE_CXX_STANDARD 17)

#******************************************************************************************  Google-tests
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)
#******************************************************************************************  Google-tests


SET(CMAKE_CXX_FLAGS "-std=c++17 -O3")
#***************************** Cluster Editing

include_directories(
        include
        include/graphs
        include/combinatorics
        include/utils
        include/datastructures
        include/clues
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

        src/clues/CETestGraphGenerator.cpp

        ) # include all source files

#***************************** Cluster Editing



add_executable(CluES "main.cpp" ${SOURCES})
target_link_libraries( CluES stdc++fs )
