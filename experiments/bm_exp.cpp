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
        return filepath1 + sep + cname + "_" + conv(psid)
                + "_" + to_string(inst_id) + ".txt";
    }

    static string filepath1;
    static streambuf* old_clog_buf;
};

string BM_CluES_exp::filepath1 = "ClusterEditingTests";

constexpr int ITERATIONS = 10;
constexpr int REPETITIONS = 1;
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

void runCluES(VVI & V, Config cnf){


}

void runOnInstance( benchmark::State & st, const string& cname ){

    int iter = 0;
    int psid = st.range(0);
    int ls_mask = st.range(1);

    for( auto _ : st ){
        st.PauseTiming();

        auto path = BM_CluES_exp::getFilePath(cname, psid, iter);
        ifstream str(path);
        auto V = readGraph(str);
        Config cnf;

        st.ResumeTiming();




        runCluES( V, cnf );

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
    auto params_set_ids = VI{ { 0, 0, 24, 0 } };
    auto delim = "__";

    constexpr bool register_ls_iteration_tests = true;

    if constexpr (register_ls_iteration_tests) {

        for (auto [cname, psid]: zip(class_names, params_set_ids)) {

            string name = cname;

            auto bm = benchmark::RegisterBenchmark(name, runOnInstance, cname)
                    ->ArgsProduct({
                        benchmark::CreateDenseRange(0, psid, 1),
                        { NM, NM+EM, NM+EM+TM, NM+CA+CR, NM+EM+CA+CR } // masks with used ls moves
                    });

            specifyBenchmarkParameters(bm);
        }
    }


//    benchmark::SetBenchmarkFilter( "A*" );
//    benchmark::SetBenchmarkFilter( "B*" );
//    benchmark::SetBenchmarkFilter( "C*" );
//    benchmark::SetBenchmarkFilter( "D*" );

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
}