//
// Created by sylwester on 2/11/24.
//

#include <filesystem>
#include <utility>
#include "StandardUtils.h"
#include "Stopwatch.h"
#include "benchmark/benchmark.h"
#include "GraphReader.h"
#include "heur/Config.h"
#include "clues/main_CE.h"
#include "heur/StateImprovers/NodeEdgeGreedyNomap.h"
#include "heur/StateImprovers/NodeEdgeGreedy.h"


class BM_CluES_exp : public benchmark::Fixture {
public:

    static void SetUpTestSuite() {
        old_clog_buf = clog.rdbuf();
        clog.rdbuf(cout.rdbuf());

    }

    static void TearDownTestSuite() {
        clog.rdbuf(old_clog_buf );
        clog << fixed;
        clog.precision(2);
    }

    static string conv( int psid ){
        stringstream str;
        if(psid < 10) str << "0";
        if( psid < 100 ) str << "0";
        str << psid;
        return str.str();
    }

    /**
     * Creates and returns a relative path for given sample name and size.
     */
    static string getFilePath(string cname, int psid, int inst_id){
        auto sep = filesystem::path::preferred_separator;
//        auto sep = "//";
        return filepath1 + sep + cname + "_" + conv(psid)
                + "_" + to_string(inst_id) + ".txt";
    }

    static string filepath1;
    static streambuf* old_clog_buf;
};

string BM_CluES_exp::filepath1 = "inputs";
//string BM_CluES_exp::filepath1 = "../experiments/ClusterEditingTests";

constexpr int ITERATIONS = 1;
constexpr int REPETITIONS = 10;
constexpr benchmark::TimeUnit unit = benchmark::TimeUnit::kMillisecond;

constexpr int NM = 1 << 0; // node move
constexpr int EM = 1 << 1; // edge move
constexpr int TM = 1 << 2; // triangle move
constexpr int CJ = 1 << 3; // cluster join
constexpr int NS = 1 << 4; // node swap
constexpr int DM = 1 << 5; // double move
constexpr int CA = 1 << 6; // component attraction
constexpr int CR = 1 << 7; // component repulsion


VVI readGraph(istream & str){
    return GraphReader::readGraphDIMACSWunweighed( str, false );
}

auto getFileComment(auto str) -> string{
    string s;
    getline(str,s);
    return s.substr( 2, s.size()-2 );
}

auto runCluESForConfiguration(VVI & V, Config cnf){

    auto counters = runCluES(V,cnf);

    return counters;
}

auto runNegForConfigurations( VVI V, Config cnf ){
    VI init_part(V.size());
    iota(ALL(init_part),0);

    ClusterGraph clg(&V,init_part);
    State st(clg, RANDOM_MATCHING);
    NEG* neg;
    constexpr bool use_neg_nomap = true;
    if constexpr (use_neg_nomap){
        neg = new NodeEdgeGreedyNomap(st);
    }else{
        neg = new NodeEdgeGreedy(st);
    }
    neg->setConfigurations(cnf);

    map<string,int64_t> counters;
    neg->perturb_mode = 0; // splitting clusters

    Global::startAlg();
    neg->counters = &counters;
    neg->improve();
    clog << "Best result found by neg: " << neg->best_result << endl;


    return counters;
}

auto createConfiguration( int rec_depth, int ls_mask, int max_perturbations,
                          int nonstandard_move_frequencies, int ls_iterations ){

    Config cnf;
    cnf.max_recursion_depth = rec_depth;
    cnf.neg_max_perturb = max_perturbations;
    if( max_perturbations == 0 ) cnf.neg_allow_perturbations = false;

    cnf.neg_max_iterations_to_do = ls_iterations;


    cnf.neg_use_edge_swaps = ls_mask & EM;
    if(ls_mask & EM){
        cnf.neg_edge_swaps_frequency = nonstandard_move_frequencies++;
    }

    cnf.neg_use_chain2_swaps = ls_mask & DM;
    if(ls_mask & DM){
        cnf.neg_chain2_swaps_frequency = nonstandard_move_frequencies++;
    }

    cnf.neg_use_triangle_swaps = ls_mask & TM;
    if(ls_mask & TM){
        cnf.neg_triangle_swaps_frequency = nonstandard_move_frequencies++;
        cnf.neg_use_triangle_swaps_to_other_clusters = false;
    }

    cnf.neg_use_join_clusters = ls_mask & CJ;
    if( ls_mask & CJ ){
        cnf.neg_join_clusters_frequency = nonstandard_move_frequencies++;
    }

    cnf.neg_use_node_interchange = ls_mask & NS;
    if( ls_mask & NS ){
        cnf.neg_node_interchanging_frequency = nonstandard_move_frequencies++;
    }



    cnf.swpCndCreatorsToUse.clear();

    cnf.use_component_repulsion = ls_mask & CA;
    if( ls_mask & CA) cnf.swpCndCreatorsToUse.push_back( exp_ord_rep );

    cnf.use_component_attraction = ls_mask & CA;
    if(ls_mask & CA) cnf.swpCndCreatorsToUse.push_back( exp_ord_attr );

    return cnf;
}

void runOnInstance( benchmark::State & st, const string& cname ){

    static map<string,int> rep_cnt;

    int psid = st.range(0);
    int ls_mask = st.range(1);
    int rec_depth = st.range(2);
    int max_perturbations = st.range(3);
    int nonstandard_move_frequencies = st.range(4);
    int ls_iterations = st.range(5);

    string rep_id = cname + "_" + to_string(psid) + "_" + to_string(ls_mask);
    int rep = rep_cnt[rep_id]++;



    for( auto _ : st ){
        st.PauseTiming();

        auto path = BM_CluES_exp::getFilePath(cname, psid, rep);
        assert( filesystem::exists(path) );// DEBUG(path);
        clog << "Running test for path: " << path << endl;
        ifstream str(path);
        auto V = readGraph(str);
        str.close();

        Config cnf = createConfiguration( rec_depth, ls_mask, max_perturbations,
                              nonstandard_move_frequencies, ls_iterations );
        st.ResumeTiming();


        map<string, int64_t> counters;
        if(rec_depth == 0){
            counters = runNegForConfigurations(V, cnf);
        }
        else{
            counters = runCluESForConfiguration( V, cnf );
        }

//        clog << "There are " << counters.size() << " counters created" << endl;

        // update counters
        st.PauseTiming();
        for( auto & [k,v] : counters ) st.counters[k] = v;
        st.ResumeTiming();

    }


}


void specifyBenchmarkParameters(auto * bm){
    bm->Unit(unit)
            ->Iterations(ITERATIONS)
            ->Repetitions(REPETITIONS)
            ->DisplayAggregatesOnly(true)
            ->UseRealTime()
            ->ComputeStatistics("max", [](const std::vector<double> &v) -> double {
                return *(std::max_element(std::begin(v), std::end(v)));
            })
            ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
                return *(std::min_element(std::begin(v), std::end(v)));
            });
}



int main(int argc, char** argv) {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    using namespace StandardUtils;


    auto class_names = vector<string>{ {"A", "B", "C", "D"} };
//    auto params_set_ids = VI{ { 22, 24, 20, 24 } };
    auto params_set_ids = VI{ { 3,3,3,3 } }; // #TEST - just to create small benchmarks for plotting
    string delim = "__";

    constexpr bool register_rec_and_perturb = false;

    if constexpr (register_rec_and_perturb) {

        for (auto [cname, psid]: zip(class_names, params_set_ids)) {

            string name = "full_clues" + delim + cname;

            auto bm = benchmark::RegisterBenchmark(name, runOnInstance, cname)
                   ->ArgsProduct({
                        benchmark::CreateDenseRange(1, psid, 1),
                        {
                            NM,
                            NM+EM,
                            NM+EM+TM,
                            NM+DM+CJ+NS,
                            NM+CA+CR,
                            NM+EM+CA+CR,
                            (1<<8) - 1 // full mask, using all moves
                        }, // masks with used ls moves
                        benchmark::CreateDenseRange(1,5,1), // recursion depth
                        {5,10,20}, // maximum number of perturbations done in NEG
                        {5,10,20}, // nonstandard move frequencies
                        {50, 100, 150}, // lS_iterations,  total number of iterations for local search
                    });
            specifyBenchmarkParameters(bm);
        }
    }





    constexpr bool register_iteration_tests = false;

    if constexpr (register_iteration_tests) {
        for (auto [cname, psid]: zip(class_names, params_set_ids)) {
            string name = "iteration_ls" + delim + cname;
            auto bm = benchmark::RegisterBenchmark(name, runOnInstance, cname)
                    ->ArgsProduct({
                                          benchmark::CreateDenseRange(1, psid, 1),
                                          {
                                                  NM,
                                                  NM + EM,
                                                  NM + EM + TM,
                                                  NM + DM + CJ + NS,
                                                  NM + CA + CR,
                                                  NM + EM + CA + CR,
                                                  (1 << 8) - 1 // full mask, using all moves
                                          }, // masks with used ls moves
                                          {0}, // recursion depth
                                          benchmark::CreateDenseRange(10, 50,
                                                                      10), // maximum number of perturbations done in NEG
                                          benchmark::CreateDenseRange(5, 50, 5), // nonstandard move frequencies
                                          {50, 100, 150}, // lS_iterations,  total number of iterations for local search
                                  });
            specifyBenchmarkParameters(bm);
        }
    }



    constexpr bool register_minimalistic_iteration_tests = true;

    if constexpr (register_minimalistic_iteration_tests) {
        for (auto [cname, psid]: zip(class_names, params_set_ids)) {
            string name = "iteration_min_ls" + delim + cname;
            auto bm = benchmark::RegisterBenchmark(name, runOnInstance, cname)
//                    ->ArgsProduct({
//                              benchmark::CreateDenseRange(1, psid, 1),
//                              {
//                                      NM,
//                                      NM+EM,
//                                      NM+EM+TM,
//                                      NM+DM+CJ+NS,
//                                      NM+EM+TM+CJ+NS+DM,
//    //                                                  (1<<8) - 1 // full mask, using all moves
//                              }, // masks with used ls moves
//                              {0}, // recursion depth
//                              {0}, // maximum number of perturbations done in NEG
//                              benchmark::CreateDenseRange(5,30,5), // nonstandard move frequencies
//                              {50, 100, 150}, // lS_iterations,  total number of iterations for local search
//                      });
                    ->ArgsProduct({
                              benchmark::CreateDenseRange(1, psid, 1),
                              {
                                      NM,
                                      NM+EM,
//                                      NM+EM+CJ+NS+DM,
//                                      NM+EM+TM,
//                                      NM+EM+TM+CJ+NS+DM,
                                      //                                                  (1<<8) - 1 // full mask, using all moves
                              }, // masks with used ls moves
                              {0}, // recursion depth
                              {0}, // maximum number of perturbations done in NEG
                              {10}, // nonstandard move frequencies
                              {100}, // ls_iterations,  total number of iterations for local search
                      });
            specifyBenchmarkParameters(bm);
        }
    }







//    benchmark::SetBenchmarkFilter( "A*" );
//    benchmark::SetBenchmarkFilter( "B*" );
//    benchmark::SetBenchmarkFilter( "C*" );
//    benchmark::SetBenchmarkFilter( "D*" );
    benchmark::SetBenchmarkFilter( "iteration_min_ls" );

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
}