
#include <clues/heur/StateImprovers/NodeEdgeGreedy.h>
#include <clues/heur/StateImprovers/NodeEdgeGreedyNomap.h>
#include <clues/heur/StateImprovers/SparseGraphTrimmer.h>
#include <clues/test_graphs.h>
#include <clues/heur/StateImprovers/NodeEdgeGreedyW1.h>
#include "clues/main_CE.h"


map<string,int64_t> runCluES( VVI V, Config cnf ){

    Global::startAlg();
    Global::counters.clear();

    Global::max_runtime_in_seconds = 300; // 5 minutes max for single test run
    clog << "Setting maximal time to " << Global::max_runtime_in_seconds << " seconds" << endl;

    TimeMeasurer::start( "Total time" );

    Global::increaseStack();

    if(!Global::disable_all_logs) GraphUtils::writeBasicGraphStatistics(V);

    VI init_part(V.size());
    iota(ALL(init_part),0);

    if(!Global::disable_all_logs){
        DEBUG(cnf.granularity_frequency);
        DEBUG(cnf.max_recursion_depth);
        DEBUG(cnf.neg_move_frequency);
        DEBUG(cnf.neg_use_edge_swaps);
        DEBUG(cnf.neg_use_triangle_swaps);
        DEBUG(cnf.neg_use_node_interchange);
        DEBUG(cnf.neg_use_queue_propagation);
        DEBUG(cnf.neg_use_join_clusters);
        DEBUG(cnf.neg_use_chain2_swaps);
        DEBUG(cnf.use_neg_map_version);
    }



    VPII best_mods; int best_result = 1e9;

    Solver solver(V, init_part, cnf);

    ClusterGraph clg(&V,init_part);
    State st(clg, RANDOM_MATCHING);
    NEG* neg = new NodeEdgeGreedyW1(st);
    neg->setConfigurations(cnf);

    neg->perturb_mode = 0; // splitting clusters
    neg->prefer_cluster_mode = 1; // prefer moving to smaller clusters
    neg->improve();

//    { // CAUTION - setting values to unused solver
//        solver.best_result = neg->best_result;
//        solver.best_partition = neg->best_partition;
//    }

    delete neg;
//    solver.run_recursive();
    solver.run_fast();

    if( solver.best_result < best_result ){
        best_result = solver.best_result;
        best_mods = solver.getModifications();
    }

    if(!Global::disable_all_logs){
        clog << "Creators: (calls,improvements):" << endl;
        for( auto & [s,p] : solver.local_search_creator_calls ){
            clog << s << " --> " << p << endl;
        }
        clog << endl << endl << endl << endl << "********************* best result found: "
                                       << best_result << endl << endl;
    }

    bool write_mods = false;
    if (write_mods) {
        VPII mods = best_mods;
        for (auto e : mods) cout << e.first+1 << " " << e.second+1 << endl;
    }

    cerr << "Final result: " << best_result << endl;
    cerr << "Total real time: " << Global::secondsFromStart() << endl;


    if(!Global::disable_all_logs) {
        TimeMeasurer::stop("Total time");
        TimeMeasurer::write();
    }

    return Global::counters;
}
