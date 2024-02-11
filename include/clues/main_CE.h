#ifndef ALGORITHMSPROJECT_MAIN_CE_H
#define ALGORITHMSPROJECT_MAIN_CE_H

#include <graphs/GraphReader.h>
#include <clues/heur/PaceUtils.h>
#include <clues/heur/ExpansionOrder.h>
#include <clues/heur/Cluster.h>
#include <clues/heur/Config.h>
#include <clues/heur/State.h>
#include <clues/heur/SwapCandidates/SwapCandidate.h>
#include <clues/heur/EOCreators/ComponentExpansion.h>
#include <clues/heur/Global.h>
#include <clues/kernelization/CEKernelizer.h>
#include <clues/kernelization/CriticalClique.h>
#include <graphs/GraphUtils.h>
#include <utils/RandomNumberGenerators.h>
#include <graphs/GraphTrimmer.h>
#include <clues/heur/Solver.h>

#include "Makros.h"

map<string,int64_t>  runCluES(VVI V, Config cnf);

#endif //ALGORITHMSPROJECT_MAIN_CE_H
