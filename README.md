# CluES - a heuristic solver for the cluster editing problem

This repository contains the source code of a heuristic cluster editing solver CluES. It also contains code used to generate instances on which CluES was compared with other state-of-the-art solvers for the problem, and supplementary material containing detailed results.

A variety of heuristics are used to find a small cluster editing set of given graph.<br>
The main algorithm works in iterations, in each iteration a new cluster editing set is found (iterations are independent from each other, the longer the solver runs, the more iterations are done and the greater chance of finding a good set).<br>


**Requirements**:

CMake VERSION 3.10.2 or higher<br>
c++ 17 or higher

<br>

**Installation**:

Use cmake to obtain a binary file, e.g. in linux in the main directory you can use the following commands:

mkdir build<br>
cd build<br>
cmake ..<br>
make

After this, the executable file named "CluES" should be in the "build" directory

<br>

**Usage:**

Given a graph in a file example_input.gr, you can run CluES in the following way
 
./CluES < example_input.gr > example_output.out 2>example_logs.err

CluES will run for exactly 600 seconds. 
Maximum run time can be set by changing variable _Global::max_runtime_in_seconds_ to a specified value in main_CE.cpp file.

<br>

**Generating tests:**

Code used to generate tests can be found on the _test-generator_ branch.

After bulding an executable file, just run it with <code> ./CluES </code> command. In the current directory a folder with tests named ClusterEditingTests will be created.