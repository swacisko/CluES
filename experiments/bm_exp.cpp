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

int PSID_FACTOR = 1;

VVI readGraph(istream & str){
    return GraphReader::readGraphDIMACSWunweighed( str, false );
}

auto getFileComment(auto str) -> string{
    string s;
    getline(str,s);
    return s.substr( 2, s.size()-2 );
}

void addStatsForCounters( auto& stats, auto &counters, auto& partition ){
    counters["res_clusters"] = set<int>(ALL(partition)).size();
    counters["res_insertions"] = get<0>(stats);
    counters["res_deletions"] = get<1>(stats);
    counters["res_modifications"] =get<2>(stats);
    counters["res_val"] = counters["res_modifications"];
}


auto runCluESForConfiguration(VVI & V, Config cnf){
    auto [counters, partition] = runCluES(V,cnf);

    auto stats = PaceUtils::getEdgeModificationStatistics( V,partition );
    addStatsForCounters(stats, counters, partition);

    map<string,int64_t> res;
    for(auto & [k,v] : counters){
        if( !k.starts_with("iter_") ) res[k]=v;
    }


    return res;
}


auto runNegForConfigurations( VVI V, Config cnf, bool use_neg_nomap = true ){
    VI init_part(V.size());
    iota(ALL(init_part),0);

    ClusterGraph clg(&V,init_part);
    State st(clg, RANDOM_MATCHING);
    NEG* neg;
    if(use_neg_nomap){
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

    auto stats = PaceUtils::getEdgeModificationStatistics( V, neg->best_partition );
    addStatsForCounters(stats, counters, neg->best_partition);

    return counters;
}

auto createConfiguration( int rec_depth, int ls_mask, int max_perturbations,
                          int nonstandard_move_frequencies, int ls_iterations, bool use_neg_nomap ){

    Config cnf;
    cnf.use_neg_map_version = use_neg_nomap;
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
    bool use_neg_nomap = st.range(6);

    string rep_id = cname + "_" + to_string(psid) + "_" + to_string(ls_mask) + to_string(rec_depth)
            + to_string(max_perturbations) + to_string(nonstandard_move_frequencies) + to_string( ls_iterations )
            + to_string(use_neg_nomap);
    int rep = rep_cnt[rep_id]++;

    psid *= PSID_FACTOR; // #TEST #CAUTION


    auto path = BM_CluES_exp::getFilePath(cname, psid, rep);
    if( !filesystem::exists(path) ) {
        st.SkipWithMessage( "File " + path + " does not exist, perhaps sizes were out of bounds"  );
        return;
    }
    clog << "Running test for path: " << path << endl;

    for( auto _ : st ){

        st.PauseTiming();

        ifstream str(path);
        auto V = readGraph(str);
        str.close();

        Config cnf = createConfiguration( rec_depth, ls_mask, max_perturbations,
                              nonstandard_move_frequencies, ls_iterations, use_neg_nomap );
        st.ResumeTiming();

        map<string, int64_t> counters;
        if(rec_depth == 0){
            counters = runNegForConfigurations(V, cnf, use_neg_nomap);
        }
        else{
            counters = runCluESForConfiguration( V, cnf );
        }

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
//    auto params_set_ids = VI{ { 22, 24, 20, 24 } }; // full tests - all instances

//    auto params_set_ids = VI{ { 9,9,9,9 } }; // #TEST - just to create small benchmarks for plotting
//    PSID_FACTOR = 2;

    auto params_set_ids = VI{ { 1,1,1,1 } }; // #TEST - t quickly check if everything works ok
    PSID_FACTOR = 1;

    string delim = "__";

    constexpr bool register_recursion_tests = true;

    if constexpr (register_recursion_tests) {

        for (auto [cname, psid]: zip(class_names, params_set_ids)) {

            string name = "recursion" + delim + cname;

            auto bm = benchmark::RegisterBenchmark(name, runOnInstance, cname)
                   ->ArgsProduct({
                        benchmark::CreateDenseRange(1, psid, 1),
                        {
                            NM,
                            NM+EM,
                            NM+CA+CR,
                            NM+EM+CA+CR,
                        }, // masks with used ls moves
                        benchmark::CreateDenseRange(1,5,1), // recursion depth, if 0 then only NEG will be used
                        {5}, // maximum number of perturbations done in NEG
                        {20}, // nonstandard move frequencies
                        {100}, // lS_iterations,  total number of iterations for local search
                        {1}// use_neg_nomap - 0 or 1
                    });
            specifyBenchmarkParameters(bm);
        }
    }

    constexpr bool register_algorithm_tests = true;

    if constexpr (register_algorithm_tests) {

        for (auto [cname, psid]: zip(class_names, params_set_ids)) {

            string name = "algorithm" + delim + cname;

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
                                      NM+EM+DM+CJ+NS+CA+CR,
                                      NM+EM+TM+DM+CJ+NS+CA+CR,
                              }, // masks with used ls moves
                              {3}, // recursion depth, if 0 then only NEG will be used
                              {5}, // maximum number of perturbations done in NEG
                              {20}, // nonstandard move frequencies
                              {100}, // ls_iterations,  total number of iterations for local search
                              {1}// use_neg_nomap - 0 or 1
                      });
            specifyBenchmarkParameters(bm);
        }
    }


    constexpr bool register_ls_iteration_tests = true;

    if constexpr (register_ls_iteration_tests) {
        for (auto [cname, psid]: zip(class_names, params_set_ids)) {
            string name = "ls_iteration" + delim + cname;
            auto bm = benchmark::RegisterBenchmark(name, runOnInstance, cname)
                    ->ArgsProduct({
                              benchmark::CreateDenseRange(1, psid, 1),
                              {
                                  NM,
                                  NM+EM,
                                  NM+EM+CJ+NS+DM,
                              }, // masks with used ls moves
                              {0}, // recursion depth, if 0 then only NEG will be used
                              {0}, // maximum number of perturbations done in NEG
                              {20}, // nonstandard move frequencies
                              {25, 50, 75, 100, 125, 150}, // ls_iterations,  total number of iterations for local search
                              {1} // use_neg_nomap
                      });
            specifyBenchmarkParameters(bm);
        }
    }


    constexpr bool register_perturbation_tests = true;

    if constexpr (register_perturbation_tests) {
        for (auto [cname, psid]: zip(class_names, params_set_ids)) {
            string name = "perturbation" + delim + cname;
            auto bm = benchmark::RegisterBenchmark(name, runOnInstance, cname)
                    ->ArgsProduct({
                              benchmark::CreateDenseRange(1, psid, 1),
                              {
                                  NM,
                                  NM+EM,
                                  NM+EM+CJ+NS+DM,
                              }, // masks with used ls moves
                              {0}, // recursion depth, if 0 then only NEG will be used
                              benchmark::CreateDenseRange(0,20,2), // maximum number of perturbations done in NEG
                              {20}, // nonstandard move frequencies
                              {(int)1e9}, // no limit for iterations, perturbations are limited
                              {0}
                      });
            specifyBenchmarkParameters(bm);
        }
    }



    constexpr bool register_neg_map_tests = true;

    if constexpr (register_neg_map_tests) {
        for (auto [cname, psid]: zip(class_names, params_set_ids)) {
            string name = "neg_map" + delim + cname;
            auto bm = benchmark::RegisterBenchmark(name, runOnInstance, cname)
                    ->ArgsProduct({
                              benchmark::CreateDenseRange(1, psid, 1),
                              {
                                  NM,
                                  NM+EM,
                                  NM+EM+CJ+NS+DM,
                              }, // masks with used ls moves
                              {0}, // recursion depth, if 0 then only NEG will be used
                              {5}, // maximum number of perturbations done in NEG
                              {20}, // nonstandard move frequencies
                              {(int)1e9}, // no limit for iterations, perturbations are limited
                              {0,1}
                      });
            specifyBenchmarkParameters(bm);
        }
    }







//    benchmark::SetBenchmarkFilter( "A*" );
//    benchmark::SetBenchmarkFilter( "B*" );
//    benchmark::SetBenchmarkFilter( "C*" );
//    benchmark::SetBenchmarkFilter( "D*" );
//    benchmark::SetBenchmarkFilter( "ls_iteartions" );
//    benchmark::SetBenchmarkFilter( "recursion" );

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
}

/**
./CluES_exp --benchmark_filter=iteration_ls__A --benchmark_fotmat=json --benchmark_out=iterls_A.json

 ./CluES_exp --benchmark_filter=perturbation__D --benchmark_fotmat=json --benchmark_out=perturbations_D.json
*/