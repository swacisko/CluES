//
// Created by sylwester on 2/6/22.
//
#include "Makros.h"
#include "CollectionOperators.h"
#include "graphs/GraphUtils.h"
#include "graphs/GraphReader.h"
#include "datastructures/FAU.h"

#include <filesystem>
namespace fs = std::filesystem;

typedef vector<string> VS;

namespace STATS {

//VS soft = { "kapoce", "clues", "mu_solver", "kanpai" }; // software to compare
    VS soft = {"kanpai", "mu_solver", "kapoce", "clues"}; // software to compare
    vector<char> test_types = {'A', 'C', 'D', 'E'};

    string test_cases_path = "../ClusterEditingTests"; // path to directory with test cases
    string soft_results_path = "results_mods"; // path to a directory with results of solvers from [soft]
    string csv_results_path = "results_full.csv";
    char sep = fs::path::preferred_separator;


/**
 * Writes strings in a separate csv line to given stream
 */
    void writeCSVLine(ostream &str, vector <string> line) {
        int cnt = 0;
        for (string s : line) {
            if (cnt++ > 0) str << ",";
            str << "\"" << s << "\"";
        }
        str << endl;
    }


    struct ClusterStats {
        int cl_cnt = 0;
        int min_cl_size = 1e9;
        int max_cl_size = -1;
        double avg_cl_size = 0;
        double cl_size_stddev = 0;

        int deletions = 0;
        int insertions = 0;
        int modifications = 0;

        // {size,cnt} - there are cnt clusters with given size. vector sorted in non-ascending order
        VPII cl_size_distr; //

        friend ostream &operator<<(ostream &str, const ClusterStats &stats);
    };

    ostream &operator<<(ostream &str, const ClusterStats &stats) {
        str << "cl_cnt: " << stats.cl_cnt << endl;
        str << "min_cl_size: " << stats.min_cl_size << endl;
        str << "max_cl_size: " << stats.max_cl_size << endl;
        str << "avg_cl_size: " << stats.avg_cl_size << endl;
        str << "cl_size_stddev: " << stats.cl_size_stddev << endl;
        str << "insertions: " << stats.insertions << endl;
        str << "deletions: " << stats.deletions << endl;
        str << "cl_size_distr: " << stats.cl_size_distr << endl;

        return str;
    }

    ClusterStats getStatsForModifications(VVI &V, VPII &mods) {
        int N = V.size();

        VPII ed = GraphUtils::getGraphEdges(V);
        for (auto &[a, b] : ed) if (a > b) swap(a, b);
        set <PII> edges(ALL(ed));

        ClusterStats stats;
        stats.modifications = mods.size();

        FAU fau(N);

        for (auto &[a, b] : mods) {
            a--;
            b--;
        }
        for (auto &[a, b] : mods) if (a > b) swap(a, b);
        set <PII> mods_zb(ALL(mods));
        for (auto &[a, b] : edges) if (mods_zb.count({a, b}) == 0) fau.Union(a, b); // not modified present edges

        for (auto &[a, b] : mods) {
            if (edges.count({a, b})) stats.deletions++;
            else {
                fau.Union(a, b);
                stats.insertions++;
            }
        }

        VVI clusters;
        {
            map<int, VI> cl_map;
            for (int i = 0; i < N; i++) cl_map[fau.Find(i)].push_back(i);
            for (auto[id, v] : cl_map) clusters.push_back(v);

            sort(ALL(clusters), [](auto a, auto b) {
                return a.size() < b.size();
            });
        }


        // now cluster stats
        assert(accumulate(ALL(clusters), 0, [](int s, VI &v) { return s + v.size(); }) == N);

        int C = clusters.size();
        stats.cl_cnt = clusters.size();
        stats.min_cl_size = clusters[0].size();
        stats.max_cl_size = clusters.back().size();
        stats.avg_cl_size = 1.0 * N / C;

        for (VI &v : clusters) stats.cl_size_stddev += (stats.avg_cl_size - v.size()) * (stats.avg_cl_size - v.size());
        stats.cl_size_stddev /= C;
        stats.cl_size_stddev = sqrt(stats.cl_size_stddev);

        map<int, int> size_cnt;
        for (auto &v : clusters) size_cnt[v.size()]++;
        for (auto[size, cnt] : size_cnt) stats.cl_size_distr.push_back({cnt, size});
        sort(ALL(stats.cl_size_distr), [](auto a, auto b) {
            if (a.first != b.first)
                return a.first > b.first;
            else return a.second > b.second;
        });

        return stats;
    }

/**
 * Finds results for each solver and writes it to given map.
 * Then   res[ {soft, test_case} ] = ClusterStats object
 */
    map <pair<string, string>, ClusterStats> getResults(VS &test_cases) {
        map <pair<string, string>, ClusterStats> res;


        int cnt = 1;
        for (auto tc : test_cases) {
            VVI V;
            { // read and create graph
                string graph_path = test_cases_path + sep + tc;
                clog << "Reading graph from file: " << graph_path << endl;
                ifstream str(graph_path);
                V = GraphReader::readGraphDIMACSWunweighed(str);
                str.close();
            }

            if (cnt++ > 3) break;

            for (auto s : soft) {
//            if( s != "clues" ){
//                clog << "Considering only CluES now" << endl;
//                continue;
//            }

                VPII mods;
                ClusterStats stats;

                { // read and create stats
                    string results_path = soft_results_path + sep + s + sep + tc;
                    clog << "Reading " << s << " result from path " << results_path << endl;
                    if (!fs::exists(results_path)) {
                        clog << "Results path " << results_path << " does not exist" << endl;
                        continue;
                    }

                    ifstream str(results_path);
                    int a, b;
                    while (str >> a) {
                        str >> b;
                        mods.push_back({a, b});
                    }
                    str.close();
                    stats = getStatsForModifications(V, mods);
                }

                res[{s, tc}] = stats;
            }
        }

        return res;
    }

/**
 * Creates statistics and writes them to a .csv file with given path [filepath]
 * This function must be called from directory when [ClusterEditingTests] directory with test exists and
 * [results] directory with results exists.
 */
    void gatherAndWriteStatistics() {

        VS test_cases;
        for (const auto &entry : fs::directory_iterator(test_cases_path)) test_cases.push_back(entry.path().filename());
        sort(ALL(test_cases));
        clog << "There are " << test_cases.size() << " test cases" << endl;

        auto results = getResults(test_cases);
        clog << "Results created" << endl;

        vector <string> present_test_cases;
        {
            set <string> zb;
            for (auto[p, stats] : results) zb.insert(p.second);
            present_test_cases = vector<string>(ALL(zb));
        }

        vector <string> present_soft;
        {
            set <string> zb;
            for (auto[p, stats] : results) zb.insert(p.first);

            for (auto s : soft) if (zb.count(s) > 0) present_soft.push_back(s);
        }


        vector <string> header = {
                "Test case",
                "N",
                "M",
                "graph avg_deg",
                "graph density",
                "test_type"
        };
        header += present_soft;


//    for( auto [p,stats] : results ){
//        clog << p << ":" << endl;
//        clog << stats << endl;
//        ENDL(1);
//    }

        VS statistics({
                              "modifications",
                              "insertions",
                              "deletions",
                              "clusters_cnt",
                              "min_cl_size",
                              "max_cl_size",
                              "avg_cl_size",
                              "cl_size_stddev",
                              "cl_size_distr"
                      });

        ofstream out_file(csv_results_path);
        out_file << fixed;
        out_file.precision(4);


        writeCSVLine(out_file, header);
        writeCSVLine(out_file, {});

        for (auto tc : present_test_cases) { // writ full statistics about modifications made by solvers
            clog << "Writing full stats for test case: " << tc << endl;

            VVI V;
            string test_type;

            { // read and create graph
                string graph_path = test_cases_path + sep + tc;
//            clog << "Reading graph from file: " << graph_path << endl;
                ifstream str(graph_path);
                V = GraphReader::readGraphDIMACSWunweighed(str);
                str.close();

                str.open(graph_path);
                getline(str, test_type);
                str.close();
            }

            int N = V.size();
            int M = GraphUtils::countEdges(V);
            double avg_deg = 2.0 * M / N;
            double dns = 2.0 * M / (1ll * N * (N - 1));

            vector <string> line = {
                    tc,
                    to_string(N),
                    to_string(M),
                    to_string(avg_deg),
                    to_string(dns),
                    test_type
            };

            for (const auto &s : present_soft) {
                auto stats = results[{s, tc}];
                line.push_back(to_string(stats.modifications));
            }

            writeCSVLine(out_file, line);
        }


        for (int i = 0; i < 5; i++) writeCSVLine(out_file, {});
        writeCSVLine(out_file, {"Full stats"});
        writeCSVLine(out_file, header);
        writeCSVLine(out_file, {});

        for (auto tc : present_test_cases) { // write full statistics about modifications made by solvers
            clog << "Writing full stats for test case: " << tc << endl;

            VVI V;
            string test_type;

            { // read and create graph
                string graph_path = test_cases_path + sep + tc;
//            clog << "Reading graph from file: " << graph_path << endl;
                ifstream str(graph_path);
                V = GraphReader::readGraphDIMACSWunweighed(str);
                str.close();

                str.open(graph_path);
                getline(str, test_type);
                str.close();
            }

            int N = V.size();
            int M = GraphUtils::countEdges(V);
            double avg_deg = 2.0 * M / N;
            double dns = 2.0 * M / (1ll * N * (N - 1));

            vector <string> line = {
                    tc,
                    to_string(N),
                    to_string(M),
                    to_string(avg_deg),
                    to_string(dns),
                    test_type
            };


            for (int opt = 0; opt < (int) statistics.size(); opt++) {
                if (opt > 0) {
                    line = {};
                    for (int i = 0; i < header.size() - present_soft.size(); i++) line.push_back("");
                }

                for (const auto &s : present_soft) {

                    auto stats = results[{s, tc}];

                    if (opt == 0) line.push_back(to_string(stats.modifications));
                    if (opt == 1) line.push_back(to_string(stats.insertions));
                    if (opt == 2) line.push_back(to_string(stats.deletions));
                    if (opt == 3) line.push_back(to_string(stats.cl_cnt));
                    if (opt == 4) line.push_back(to_string(stats.min_cl_size));
                    if (opt == 5) line.push_back(to_string(stats.max_cl_size));
                    if (opt == 6) line.push_back(to_string(stats.avg_cl_size));
                    if (opt == 7) line.push_back(to_string(stats.cl_size_stddev));
                    if (opt == 8) {
                        stringstream str;
                        str << fixed;
                        str.precision(4);
                        for (auto p : stats.cl_size_distr) str << p << "  ";
                        line.push_back(str.str());
                    }
                }

                line.push_back("");
                line.push_back(statistics[opt]);
                writeCSVLine(out_file, line);
            }

            writeCSVLine(out_file, {});
        }


        writeCSVLine(out_file, {});
        writeCSVLine(out_file, {});


        writeCSVLine(out_file, {"st | nd | rd | th"});

        header = {"test type"};
        header += present_soft;
        writeCSVLine(out_file, header);

        for (char type : test_types) { // creating rankings on graph classes - how many times each solver took a place
            string tp = "0";
            tp[0] = type;
            vector <string> line = {tp};
            map <string, VI> places;
            for (auto s : present_soft) places[s] = VI(4, 0);

            for (auto tc : present_test_cases) {
                if (tc[0] != type) continue;

                vector <pair<string, int>> res;
                for (auto s : present_soft) res.emplace_back(s, results[{s, tc}].modifications);

                sort(ALL(res), [](auto a, auto b) {
                    return a.second < b.second;
                });

                int pl = 0;
                for (int i = 0; i < res.size(); i++) {
                    if (i > 0 && res[i].second > res[i - 1].second) pl++;
                    places[res[i].first][pl]++;
                }
            }

            for (auto s : present_soft) {
                string temp;
                for (int pl = 0; pl < 4; pl++) {
                    if (pl > 0) temp += " | ";
                    temp += to_string(places[s][pl]);
                }

                line.push_back(temp);
            }

            writeCSVLine(out_file, line);
        }


        writeCSVLine(out_file, {});
        writeCSVLine(out_file, {});
        writeCSVLine(out_file, {"PACE score"});

        header = {"test type"};
        header += present_soft;
        writeCSVLine(out_file, header);

        map<string, double> total_score;
        map<char, int> tc_count;

        for (char type : test_types) { // creating scores for graph classes
            string tp = "0";
            tp[0] = type;
            vector <string> line = {tp};

            map<pair<char, string>, double> scores;

            for (auto tc : present_test_cases) {
                if (tc[0] != type) continue;
                tc_count[type]++;

                int m = 1e9;
                for (auto s : present_soft) {
                    auto stats = results[{s, tc}];
                    m = min(m, stats.modifications);
                }

                for (auto s : present_soft) {
                    auto stats = results[{s, tc}];
                    scores[{type, s}] += 1.0 * m / stats.modifications;
                    total_score[s] += 1.0 * m / stats.modifications;
                }
            }

            for (auto s : present_soft) {
                double score = scores[{type, s}] * 100.0 / tc_count[type];
                stringstream str;
                str << fixed;
                str.precision(4);
                str << score;
                line.push_back(str.str());
            }

            writeCSVLine(out_file, line);
        }

        {// writing total score
            vector <string> line = {"Total score"};
            int total_cases = present_test_cases.size();
            for (auto s : present_soft) {
                double score = total_score[s] * 100.0 / total_cases;
                stringstream str;
                str << fixed;
                str.precision(4);
                str << score;
                line.push_back(str.str());
            }
            writeCSVLine(out_file, line);
        }


        out_file.close();

    }
};

int main_gather_stats(){

    STATS::gatherAndWriteStatistics();

    return 0;
}