//
// Created by sylwester on 1/17/22.
//

#include <graphs/generators/GraphGenerator.h>
#include <graphs/GraphUtils.h>
#include <utils/RandomNumberGenerators.h>
#include <filesystem>
#include "clues/CETestGraphGenerator.h"
#include "StandardUtils.h"
#include "CombinatoricUtils.h"
#include <filesystem>

void CETestGraphGenerator::generateAllInstances4() {
    if(out_directory == ""){
        clog << "Output directory empty, not creating instances" << endl;
        return;
    }

    using namespace std::filesystem;

    if( !is_directory(out_directory) ){
        auto pth = path(out_directory);
        create_directory(pth);
    }else{
        clog << "Directory " << out_directory << " already exists, rewriting all files" << endl;
    }

    const int MAX_E = 3'000'000;
    const int MIN_E = 10'000;

    char sep = (char)std::filesystem::path::preferred_separator;

    auto getTestName = [](int id, int instance){
        if( id < 10 ) return "00" + to_string(id) + "_" + to_string(instance);
        else if( id < 100 ) return "0" + to_string(id) + "_" + to_string(instance);
        else return to_string(id) + "_" + to_string(instance);;
    };


    const bool generate_nonuniform = true;
    const bool generate_clustered_graphs = true;
    const bool generate_spatial_graphs = true;
    const bool generate_kpartite_graphs = true;

    const int INSTANCES_PER_GRAPH_CLASS = 10;

    if(generate_nonuniform){ // nonuniform graphs
        int test_id = 1;
        clog << "Generating nonuniform random graphs" << endl;

        VI nodes = {1'000, 2'000, 3'000, 4'000, 5'000};
        VD densities = {0.05, 0.1, 0.2, 0.3, 0.4};

        for(auto N : nodes) {
            for( auto dns : densities ) {

                LL M = (1ll * N * (N-1) ) >> 1;
                M *= dns;
                if( M > MAX_E ) continue;

                for( int t = 0; t < INSTANCES_PER_GRAPH_CLASS; t++ ) {

                    clog << "\rCreating test #" << test_id << flush;
                    ofstream str(out_directory + sep + "A_" + getTestName(test_id, t) + ".txt");

                    VI probabs(N);
                    iota(ALL(probabs), sqrt(N));
                    createRandomNonuniform(str, N, M, probabs);
                    str.close();
                }

                test_id++;
            }
        }

        clog << endl << "Nonuniform tests generated" << endl;
        ENDL(3);
    }

    //********************************************************************************************

    /**
     * When creating cluster sizes or partition sizes in some way (e.g. randomly), it may be impossible to
     * create a graph with exactly predefined number of nodes N. Thus, we at the end smooth sizes - we increase
     * smallest clusters / decrease largest clusters until the graph size is exactly N.
     *
     * #CAUTION! It may be poossible (but with exponentially decreasing probability), that it will not be possible
     * to smooth the graph with current implementation and the algorithm will get stuck in infinite loop.
     */
    auto smoothClusterSizes = []( VI& cl_sizes, int N, int minC, int maxC ){
        // creating cluster sizes and 'smoothing' them to exactly given N
        int K = cl_sizes.size();
        int totN = accumulate(ALL(cl_sizes),0);

        int p = 0;
        int AVG_S = (minC + maxC) / 2;
        int cnt = 0;
        bool opt = false;

        while( totN > N ){
            if( opt || cl_sizes[p] > AVG_S ){
                cl_sizes[p]--;
                totN--;
            }
            p = (p+1)%K;

            if(cnt++ > 2*N) opt = true;
        }

        opt = false;
        cnt = 0;

        while( totN < N ){
            if( opt || cl_sizes[p] < AVG_S ){
                cl_sizes[p]++;
                totN++;
            }
            p = (p+1)%K;

            if(cnt++ > 2*N) opt = true;
        }
    };

    if(generate_clustered_graphs){ // C
        int test_id = 1;
        clog << "Generating clustered tests" << endl;

        VI nodes = {1'000, 2'000, 3'000};
        VI ks = {30, 100};
        vector<pair<double,double>> psqs = {
                {0.6, 0.5},
                {0.6, 0.4},
                {0.5, 0.4},
                {0.5, 0.3}
        };

        int minC, maxC;

        for(auto N : nodes){
            for( int K : ks ) {
                for (auto [p,q] : psqs) {

                    int avg_cl_size = N / K;
                    minC = avg_cl_size * 0.75;
                    maxC = avg_cl_size * 1.25;

                    for(  int t=0; t<INSTANCES_PER_GRAPH_CLASS; t++) {
                        VI cl_sizes;

                        UniformIntGenerator rnd(minC, maxC);
                        for(int i=0; i<K; i++) cl_sizes.push_back( rnd.rand() );
                        smoothClusterSizes(cl_sizes, N, minC, maxC);
                        sort(ALL(cl_sizes));

                        assert( accumulate(ALL(cl_sizes),0) == N );

                        int M = 0;
                        int UB = 0;

                        for (int i = 0; i < K; i++) {
                            int C1 = cl_sizes[i];
                            for (int j = i; j < K; j++) {
                                int C2 = cl_sizes[j];

                                if (i == j){
                                    int e = round( p * int((C1 * (C1 - 1) / 2)) );
                                    M += e;
                                    UB += (C1 * (C1-1) / 2) - e;
                                }
                                else {
                                    int e = round( q * (C1*C2) );
                                    M += e;
                                    UB += e;
                                }
                            }
                        }

                        clog << "\rCreating test #" << test_id << flush;

                        ofstream str(out_directory + sep + "C_" + getTestName(test_id, t) + ".txt");
                        createRandomClustered1(str, cl_sizes, p, q, UB, M);
                        str.close();
                    }

                    test_id++;

                }
            }
        }

        clog << endl << "clustered tests generated" << endl;
        ENDL(3);
    }



    //********************************************************************************************


    if(generate_spatial_graphs){ // D
        int test_id = 1;
        clog << "Generating spatial tests" << endl;

        VI nodes = {5'000, 10'000, 15'000, 20'000};
        VI ks = { 20,40,60,80,100 };

        for( int N : nodes ){
            for( auto K : ks ){

                LL Mmax = (N*K);
                int rand_neigh=K;

                for( int t=0; t<INSTANCES_PER_GRAPH_CLASS; t++ ) {
                    clog << "\rCreating test #" << test_id << endl;
                    ofstream str(out_directory + sep + "D_" + getTestName(test_id, t) + ".txt");
                    createRandomSpatialEuclidTorus(str, N, K, rand_neigh, 1);
                    str.close();
                }

                test_id++;
            }
        }

        clog << endl << "spatial tests generated" << endl;
        ENDL(3);
    }


    //********************************************************************************************


    if(generate_kpartite_graphs){ // E
        int test_id = 1;
        clog << "Generating k-partite tests" << endl;

        VI nodes = {1'000, 2'000, 3'000};
        VI ks = {30, 100};
        VD densities = {0.1, 0.2, 0.3, 0.4};

        for( auto N : nodes ){
            for(auto k : ks) {
                for( auto dns : densities ){
                    for( int t=0; t<INSTANCES_PER_GRAPH_CLASS; t++ ) {
                        VI cl_sizes;
                        int AVG_S = N / k;
                        int minC = 0.75 * AVG_S;
                        int maxC = 1.25 * AVG_S;

                        UniformIntGenerator rnd(minC, maxC);
                        for(int i=0; i<k; i++) cl_sizes.push_back( rnd.rand() );
                        smoothClusterSizes(cl_sizes, N, minC, maxC);
                        sort(ALL(cl_sizes));

                        assert( N == accumulate(ALL(cl_sizes),0) );
                        assert(cl_sizes.size() == k);

                        LL M = 0;
                        for( int i=0; i<k; i++ ){
                            for( int j=i+1; j<k; j++ ){
                                M += cl_sizes[i] * cl_sizes[j];
                            }
                        }

                        M *= dns;

                        if (M > MAX_E) continue;

                        clog << "\rCreating test #" << test_id << flush;
                        ofstream str(out_directory + sep + "E_" + getTestName(test_id, t) + ".txt");
                        createRandomKPartite(str, M, dns, cl_sizes);
                        str.close();
                    }

                    test_id++;
                }

            }
        }

        clog << endl << "k-partite tests generated" << endl;
        ENDL(3);
    }


    clog << "Tests generated!" << endl;
}


void CETestGraphGenerator::generateAllInstances3() {
    if(out_directory == ""){
        clog << "Output directory empty, not creating instances" << endl;
        return;
    }

    using namespace std::filesystem;

    if( !is_directory(out_directory) ){
        auto pth = path(out_directory);
        create_directory(pth);
    }else{
        clog << "Directory " << out_directory << " already exists, rewriting all files" << endl;
    }

    const int MAX_E = 3'000'000;
    const int MIN_E = 10'000;

    char sep = (char)std::filesystem::path::preferred_separator;

    auto getTestName = [](int id){
        if( id < 10 ) return "00" + to_string(id);
        else if( id < 100 ) return "0" + to_string(id);
        else return to_string(id);
    };


    const bool generate_ABF = false;
    const bool generate_clustered_graphs = false;
    const bool generate_spatial_graphs = false;
    const bool generate_kpartite_graphs = false;

    if(generate_ABF){ // ABF
        int test_id = 1;
        clog << "Generating createRandom, regular and nonuniform tests" << endl;

        VI nodes = {
                1'000,
                2'000,
                3'000,
                4'000,
                5'000
        };

        VD densities = {
                0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.45
        };

        int opt = 0;
        const int RANDOM = 0;
        const int REGULAR = 1;
        const int NONUNIFORM = 2;

        for(int opt = 0; opt < 3; opt++) {
            for(auto N : nodes) {
                for( auto dns : densities ) {

                    LL M = (1ll * N * (N-1) ) >> 1;
                    M *= dns;
    //                if( M > MAX_E || M < MIN_E ) continue;
                    if( M > MAX_E ) continue;

                    switch (opt) {
                        case RANDOM: {
                            clog << "\rCreating test #" << test_id << flush;
                            //                        ofstream str(out_directory + sep + "AA_" + getTestName(test_id) + ".txt");
                            ofstream str(out_directory + sep + "A_" + getTestName(test_id) + ".txt");
                            createRandom2(str, N, M);
                            str.close();

                            test_id++;
                            break;
                        }

                        case REGULAR: {
                            clog << "\rCreating test #" << test_id << flush;
                            //                        ofstream str(out_directory + sep + "AB_" + getTestName(test_id) + ".txt");
                            ofstream str(out_directory + sep + "A_" + getTestName(test_id) + ".txt");
                            int K = 2ll * M / N;
                            createRandomRegular(str, N, K);
                            str.close();

                            test_id++;
                            break;
                        }

                        case NONUNIFORM: {
                            if (dns <= 0.4) { // creating for dense graphs may be extremely time-consuming

                                clog << "\rCreating test #" << test_id << flush;
                                //                            ofstream str(out_directory + sep + "AF_" + getTestName(test_id) + ".txt");
                                ofstream str(out_directory + sep + "A_" + getTestName(test_id) + ".txt");

                                VI probabs(N);
                                iota(ALL(probabs), sqrt(N));
                                createRandomNonuniform(str, N, M, probabs);
                                str.close();

                                test_id++;
                                break;
                            }
                        }
                    }
                }

            }
        }

        clog << endl << "ABF tests generated" << endl;
        ENDL(3);
    }

    //********************************************************************************************

    /**
     * When creating cluster sizes or partition sizes in some way (e.g. randomly), it may be impossible to
     * create a graph with exactly predefined number of nodes N. Thus, we at the end smooth sizes - we increase
     * smallest clusters / decrease largest clusters until the graph size is exactly N.
     *
     * #CAUTION! It may be poossible (but with exponentially decreasing probability), that it will not be possible
     * to smooth the graph with current implementation and the algorithm will get stuck in infinite loop.
     */
    auto smoothClusterSizes = []( VI& cl_sizes, int N, int minC, int maxC ){
//        clog << "Smoothing clusters" << endl;
     // creating cluster sizes and 'smoothing' them to exactly given N
        int K = cl_sizes.size();
        int totN = accumulate(ALL(cl_sizes),0);

        int p = 0;
        int AVG_S = (minC + maxC) / 2;
        int cnt = 0;
        bool opt = false;

        while( totN > N ){
            if( opt || cl_sizes[p] > AVG_S ){
                cl_sizes[p]--;
                totN--;
            }
            p = (p+1)%K;

            if(cnt++ > 2*N) opt = true;
        }

        opt = false;
        cnt = 0;

        while( totN < N ){
            if( opt || cl_sizes[p] < AVG_S ){
                cl_sizes[p]++;
                totN++;
            }
            p = (p+1)%K;

            if(cnt++ > 2*N) opt = true;
        }

//        clog << "clusters smoothed" << endl;
    };

    if(generate_clustered_graphs){ // C
        int test_id = 1;
        clog << "Generating clustered tests" << endl;

        VI nodes = {3'000};
        VI ks = {30};
        vector<pair<double,double>> psqs = {
                {0.55, 0.45},
                {0.55, 0.35},
                {0.45, 0.35},
                {0.45, 0.25},
                {0.35, 0.25},
                {0.35, 0.15}
        };

        const int minC = 80, maxC = 120;

        for(auto N : nodes){
            for( int K : ks ) {
                for (auto [p,q] : psqs) {

                    for(int T=0; T<10; T++) {
                        VI cl_sizes;

                        UniformIntGenerator rnd(minC, maxC);
                        for(int i=0; i<K; i++) cl_sizes.push_back( rnd.rand() );
                        smoothClusterSizes(cl_sizes, N, minC, maxC);
                        sort(ALL(cl_sizes));


                        assert( accumulate(ALL(cl_sizes),0) == N );

                        int M = 0;
                        int UB = 0;

                        for (int i = 0; i < K; i++) {
                            int C1 = cl_sizes[i];
                            for (int j = i; j < K; j++) {
                                int C2 = cl_sizes[j];

                                if (i == j){
                                    int e = round( p * int((C1 * (C1 - 1) / 2)) );
                                    M += e;
                                    UB += (C1 * (C1-1) / 2) - e;
                                }
                                else {
                                    int e = round( q * (C1*C2) );
                                    M += e;
                                    UB += e;
                                }
                            }
                        }

//                        if (M < MIN_E || M > MAX_E) continue;

                        clog << "\rCreating test #" << test_id << flush;

                        ofstream str(out_directory + sep + "C_" + getTestName(test_id) + ".txt");
                        createRandomClustered1(str, cl_sizes, p, q, UB, M);
                        str.close();
                        test_id++;

                    }

                }
            }
        }

        clog << endl << "clustered tests generated" << endl;
        ENDL(3);
    }



    //********************************************************************************************


    if(generate_spatial_graphs){ // D
        int test_id = 1;
        clog << "Generating spatial tests" << endl;

        VI nodes = {5'000, 10'000, 15'000, 20'000};
        VI ks = { 20,40,60,80,100 };

        for( int N : nodes ){
            for( auto K : ks ){

                LL Mmax = (N*K);

//                if( Mmax > MAX_E ) continue;
//                if( Mmax/2 < MIN_E ) continue;

                int rand_neigh=K;
                clog << "\rCreating test #" << test_id << endl;
                ofstream str(out_directory + sep + "D_" + getTestName(test_id) + ".txt");
                createRandomSpatialEuclidTorus(str, N, K, rand_neigh, 1);
                str.close();
                test_id++;
            }
        }

        clog << endl << "spatial tests generated" << endl;
        ENDL(3);
    }


    //********************************************************************************************


    if(generate_kpartite_graphs){ // E
        int test_id = 1;
        clog << "Generating k-partite tests" << endl;

        VI nodes = {1'000, 2'000, 3'000, 4'000, 5'000};
        VI ks = {3, 5, 10, 20, 50, 75, 100, 150};
        VD densities = {0.1, 0.2, 0.3};

        for( auto N : nodes ){
            for(auto k : ks) {
                for( auto dns : densities ){

                    VI cl_sizes;
                    int AVG_S = N / k;
                    int minC = 0.8 * AVG_S;
                    int maxC = 1.2 * AVG_S;

                    UniformIntGenerator rnd(minC, maxC);
                    for(int i=0; i<k; i++) cl_sizes.push_back( rnd.rand() );
                    smoothClusterSizes(cl_sizes, N, minC, maxC);
                    sort(ALL(cl_sizes));

                    assert( N == accumulate(ALL(cl_sizes),0) );
                    assert(cl_sizes.size() == k);

                    LL M = 0;
                    for( int i=0; i<k; i++ ){
                        for( int j=i+1; j<k; j++ ){
                            M += cl_sizes[i] * cl_sizes[j];
                        }
                    }

                    M *= dns;

//                    if (M < MIN_E) continue;
                    if (M > MAX_E) continue;


                    clog << "\rCreating test #" << test_id << flush;
                    ofstream str(out_directory + sep + "E_" + getTestName(test_id) + ".txt");
                    createRandomKPartite(str, M, dns, cl_sizes);

                    str.close();
                    test_id++;
                }

            }
        }

        clog << endl << "k-partite tests generated" << endl;
        ENDL(3);
    }


    clog << "Tests generated!" << endl;
}






void CETestGraphGenerator::generateAllInstances2() {
    if(out_directory == ""){
        clog << "Output directory empty, not creating instances" << endl;
        return;
    }

    using namespace std::filesystem;

    if( !is_directory(out_directory) ){
        auto pth = path(out_directory);
        create_directory(pth);
    }else{
        clog << "Directory " << out_directory << " already exists, rewriting all files" << endl;
    }

    const int MAX_E = 3'000'000;
    const int MIN_E = 10'000;

    char sep = (char)std::filesystem::path::preferred_separator;

    auto getTestName = [](int id){
        if( id < 10 ) return "00" + to_string(id);
        else if( id < 100 ) return "0" + to_string(id);
        else return to_string(id);
    };


    const bool generate_ABF = false;
    const bool generate_clustered_graphs = true;
    const bool generate_spatial_graphs = false;
    const bool generate_kpartite_graphs = false;

    if(generate_ABF){ // ABF
        int test_id = 1;
        clog << "Generating createRandom, regular and nonuniform tests" << endl;

        VI nodes = {
                1'000,
                2'000,
                5'000,
                10'000
        };

        VD densities = {
            0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.47, 0.5, 0.53
        };

        int opt = 0;
        const int RANDOM = 0;
        const int REGULAR = 1;
        const int NONUNIFORM = 2;

        for(auto N : nodes) {
            for( auto dns : densities ) {

                LL M = (1ll * N * (N-1) ) >> 1;
                M *= dns;
                if( M > MAX_E || M < MIN_E ) continue;

                for(int opt = 0; opt < 3; opt++) {
                    switch (opt) {
                        case RANDOM: {
                            clog << "\rCreating test #" << test_id << flush;
                            //                        ofstream str(out_directory + sep + "AA_" + getTestName(test_id) + ".txt");
                            ofstream str(out_directory + sep + "A_" + getTestName(test_id) + ".txt");
                            createRandom2(str, N, M);
                            str.close();

                            test_id++;
                            break;
                        }

                        case REGULAR: {
                            clog << "\rCreating test #" << test_id << flush;
    //                        ofstream str(out_directory + sep + "AB_" + getTestName(test_id) + ".txt");
                            ofstream str(out_directory + sep + "A_" + getTestName(test_id) + ".txt");
                            int K = 2ll * M / N;
                            createRandomRegular(str, N, K);
                            str.close();

                            test_id++;
                            break;
                        }

                        case NONUNIFORM: {
                            if (dns <= 0.4) { // creating for dense graphs may be extremely time-consuming

                                clog << "\rCreating test #" << test_id << flush;
    //                            ofstream str(out_directory + sep + "AF_" + getTestName(test_id) + ".txt");
                                ofstream str(out_directory + sep + "A_" + getTestName(test_id) + ".txt");

                                VI probabs(N);
                                iota(ALL(probabs), sqrt(N));
                                createRandomNonuniform(str, N, M, probabs);
                                str.close();

                                test_id++;
                                break;
                            }
                        }
                    }
                }

                opt = (opt+1)%3;
            }
        }

        clog << endl << "ABF tests generated" << endl;
        ENDL(3);
    }

    //********************************************************************************************


    if(generate_clustered_graphs){ // C
        int test_id = 1;
        clog << "Generating clustered tests" << endl;

        VI nodes = {
                1000,
                2'000,
                3'000
        };

        VI ks = {
                10, 30, 100
        };

        VD qs = {
                0.2, 0.325, 0.45
        };

        VD ps = {
                0.55, 0.65, 0.75
        };

        bool opt = false;

        for(auto N : nodes){
            for( int K : ks ) {
                for (double p : ps) {

                    for (auto q : qs) {
                        VI cl_sizes;

                        if (opt) {
                            int c = (2*N/K - (K-1)) / 2;
                            for(int i=0; i<K; i++){
                                cl_sizes.push_back(c);
                                c++;
                            }
                        } else cl_sizes = VI(K, N / K);

                        if( cl_sizes[0] < 5 ) cl_sizes = VI(K, N / K);

                        int M = 0;
                        for( int i=0; i<K; i++ ){
                            int C1 = cl_sizes[i];
                            for( int j=i; j<K; j++ ){
                                int C2 = cl_sizes[j];

                                if( i == j ) M += p * (C1 * (C1-1) / 2);
                                else{
                                    M += q * ( C1 * C2 / 2 );
                                }
                            }
                        }

                        if( M < MIN_E || M > MAX_E ) continue;

                        clog << "\rCreating test #" << test_id << flush;

//                        DEBUG(cl_sizes);

                        ofstream str(out_directory + sep + "C_" + getTestName(test_id) + ".txt");
                        createRandomClustered1(str, cl_sizes, p, q, -1, -1);
                        str.close();
                        test_id++;

                        opt = !opt;
                    }
                }
            }
        }

        clog << endl << "clustered tests generated" << endl;
        ENDL(3);
    }



    //********************************************************************************************


    if(generate_spatial_graphs){ // D
        int test_id = 1;
        clog << "Generating spatial tests" << endl;

        VI nodes = {
                500,
                1000,
                2'000,
                5'000,
                10'000,
                20'000
        };

        VD densities = {
                0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7
        };

        for( int N : nodes ){
            for( auto dns : densities ){
                LL K = N * dns;

                LL Mmax = (N*K);

                if( Mmax > MAX_E ) continue;
                if( Mmax/2 < MIN_E ) continue;

                int rand_neigh=K;
                clog << "\rCreating test #" << test_id << endl;
                ofstream str(out_directory + sep + "D_" + getTestName(test_id) + ".txt");
                createRandomSpatialEuclidTorus(str, N, K, rand_neigh, 1);
                str.close();
                test_id++;
            }
        }

        clog << endl << "spatial tests generated" << endl;
        ENDL(3);
    }


    //********************************************************************************************


    if(generate_kpartite_graphs){ // E
        int test_id = 1;
        clog << "Generating k-partite tests" << endl;

        VI nodes = {
                1'000,
                2'000,
                5'000,
                10'000
        };

        VD factors = {
                1.0, 1.5
        };

        VI ks = {
                2,3,5,10
        };

        VD densities = {
//                0.01, 0.02, 0.05, 0.1,
                0.2, 0.3, 0.4, 0.5, 0.6, 0.7
        };

        int opt = 0;

        for( auto N : nodes ){
            for( auto dns : densities ){
                for(auto k : ks) {

                    auto f = factors[(opt++)&1];
                    VI cl_sizes;
                    if(f > 1.0) {
                        int n1 = N * (1.0 - f) / (1.0 - pow(f, k));
                        for (int i = 0; i < k; i++) {
                            cl_sizes.push_back(n1);
                            n1 *= f;
                        }
                    }else cl_sizes = VI( k, N/k );

                    assert(cl_sizes.size() == k);

                    LL M = 0;
                    for( int i=0; i<k; i++ ){
                        for( int j=i+1; j<k; j++ ){
                            M += cl_sizes[i] * cl_sizes[j];
                        }
                    }

                    M *= dns;

                    if (M < MIN_E) continue;
                    if (M > MAX_E) continue;


                    clog << "\rCreating test #" << test_id << flush;
                    ofstream str(out_directory + sep + "E_" + getTestName(test_id) + ".txt");
                    createRandomKPartite(str, M, dns, cl_sizes);

                    str.close();
                    test_id++;
                }

            }
        }

        clog << endl << "k-partite tests generated" << endl;
        ENDL(3);
    }


    clog << "Tests generated!" << endl;
}

void CETestGraphGenerator::generateAllInstances() {
    if(out_directory == ""){
        clog << "Output directory empty, not creating instances" << endl;
        return;
    }

    using namespace std::filesystem;

    if( !is_directory(out_directory) ){
        auto pth = path(out_directory);
        create_directory(pth);
    }else{
        clog << "Directory " << out_directory << " already exists, rewriting all files" << endl;
    }

    const int MAX_E = 2'000'000;

    char sep = (char)std::filesystem::path::preferred_separator;

    auto getTestName = [](int id){
        if( id < 10 ) return "00" + to_string(id);
        else if( id < 100 ) return "0" + to_string(id);
        else return to_string(id);
    };


    const bool generate_random_graphs = true;
    const bool generate_random_regular_graphs = true;
    const bool generate_clustered_graphs = true;
    const bool generate_spatial_graphs = true;
    const bool generate_kpartite_graphs = true;
    const bool generate_nonuniform_graphs = true;

    if(generate_random_graphs){ // A
        int test_id = 1;
        clog << "Generating createRandom2 tests" << endl;

        VI nodes = {
//                500,
                2'000, 5'000, 10'000, 20'000, 50'000
        };

        VI degrees = {
                5,10,20,50,100,200
        };


        VPII random_sizes = StandardUtils::productIf(nodes, degrees, [&](int a, int b){
            return a*b <= MAX_E;
        });

        for( auto & [N,M] : random_sizes ) M *= N;

        for( int n : VI({ 1'000, 2'000 }) ){
            LL e = 1ll*n*(n-1)/4;

            if( 1.05 * e > MAX_E ) continue;

            random_sizes.emplace_back( n, e*0.99 ); // adding fractionally less than half edges
//            random_sizes.emplace_back( n, e*1.01 ); // adding fractionally more than half edges

            random_sizes.emplace_back( n, e*0.97 ); // adding fractionally less than half edges
//            random_sizes.emplace_back( n, e*1.03 ); // adding fractionally more than half edges
        }

        for(auto[N, M] : random_sizes) {
            clog << "\rCreating test #" << test_id << flush;
            ofstream str(out_directory + sep + "A_" + getTestName(test_id) + ".txt");
            createRandom2(str, N, M);
            str.close();
            test_id++;
        }

        clog << endl << "createRandom2 tests generated" << endl;
        ENDL(3);
    }


    //********************************************************************************************


    if(generate_random_regular_graphs){ // B
        int test_id = 1;
        clog << "Generating regular tests" << endl;

        VI nodes = {
//                500,
                2'000, 5'000, 10'000, 20'000, 50'000
        };

        VI degrees = {
                5,10,20,50,100,200
        };


        VPII regular_sizes = StandardUtils::productIf(nodes, degrees, [&](int a, int b){
            return b < a && a*b/2 <= MAX_E;
        });

        for( int n : VI({1000, 2'000}) ){
            LL e = 1ll*n*(n-1)/4;
            if( 1.05*e > MAX_E ) continue;

            regular_sizes.emplace_back( n, 0.97 * n/2 );
            regular_sizes.emplace_back( n, 0.99 * n/2 );
        }

        for (auto [N, K] : regular_sizes) {
            clog << "\rCreating test #" << test_id << flush;

            ofstream str(out_directory + sep + "B_" + getTestName(test_id) + ".txt");
            createRandomRegular(str, N, K);
            str.close();
            test_id++;
        }

        clog << endl << "regular tests generated" << endl;
        ENDL(3);
    }

    //********************************************************************************************


    if(generate_clustered_graphs){ // C
        int test_id = 1;
        clog << "Generating clustered tests" << endl;

        VI clusters_cnt = {
            10,20,50,100,200
        };

        VPII cluster_size_bound = {
                {1,10},
                {10, 10},
                {10, 50},
                {20, 20},
                {50, 50},
                {1,100}
        };

        vector<pair<double,double>> probabs = {
                {0.51, 0.49},
                {0.53, 0.47}
//                {0.55, 0.45}
        };

        for( int K : clusters_cnt ){
            for( auto [minC, maxC] : cluster_size_bound ){
                for( auto [p,q] : probabs ){
                    LL avg_c_size = (minC+maxC) / 2;
                    LL expN = K * avg_c_size;

                    if( expN < 500 ) continue;

                    LL expM = K * p * avg_c_size * ( avg_c_size - 1 ) / 2
                            +  q * expN * ( expN - avg_c_size ) / 2;
                    if( 1.05 * expM > MAX_E ) continue;

                    clog << "\rCreating test #" << test_id << flush;

                    ofstream str(out_directory + sep + "C_" + getTestName(test_id) + ".txt");
                    createRandomClustered1(str, K, p, q, minC, maxC);
                    str.close();
                    test_id++;
                }
            }
        }

        clog << endl << "clustered tests generated" << endl;
        ENDL(3);
    }



    //********************************************************************************************


    if(generate_spatial_graphs){ // D
        int test_id = 1;
        clog << "Generating spatial tests" << endl;

        VI nodes = {
//                5'000, 20'000, 50'000
                10'000, 20'000, 30'000
        };

        VI ks = {
//                10,20,50,100,200
                5,10,20,50,100,200
        };

        VD sparsities = {
                1.0, 0.5
        };

        for( int N : nodes ){
            for( int K : ks ){
                for( double sparsity : sparsities ){
//                    {
//                        int expE = N * K * 0.75;
//                        if(expE > MAX_E) continue;
//                    }

                    int rand_neigh= K;

                    /*if(N == 30'000) { // usual euclidean metric only for case N = 20'000
                        clog << "\rCreating test #" << test_id << endl;
                        ofstream str(out_directory + sep + "D_" + getTestName(test_id) + ".txt");
                        createRandomSpatialEuclid(str, N, K, rand_neigh, sparsity);
                        str.close();
                        test_id++;
                    }*/

                    { // euclidean on torus for all test cases
                        clog << "\rCreating test #" << test_id << endl;
                        ofstream str(out_directory + sep + "D_" + getTestName(test_id) + ".txt");
                        createRandomSpatialEuclidTorus(str, N, K, rand_neigh, sparsity);
                        str.close();
                        test_id++;
                    }
                }
            }
        }

        clog << endl << "spatial tests generated" << endl;
        ENDL(3);
    }


    //********************************************************************************************


    if(generate_kpartite_graphs){ // E
        int test_id = 1;
        clog << "Generating k-partite tests" << endl;

        VI ks = {
                20,50,100,200,300
        };

        VPII bounds = {
                {3,10},
                {20, 20},
                {10, 100},
                {50,100}
        };

        VD probabs = {
                0.6, 0.5, 0.4
        };

        for( auto [minC, maxC] : bounds ){
            for( auto p : probabs ){
                for(auto K : ks){
                    LL avg_c_size = ((minC+maxC)>>1);
                    LL expN = K*avg_c_size;


                    LL M = expN * (expN - avg_c_size) / 2;
                    M *= p;

                    if( M < 3'000 ) continue;
                    if( M > MAX_E) continue;

                    clog << "\rCreating test #" << test_id << flush;
                    ofstream str(out_directory + sep + "E_" + getTestName(test_id) + ".txt");
                    createRandomKPartite(str, M, K,p, minC, maxC);
                    str.close();
                    test_id++;
                }

            }
        }

        clog << endl << "k-partite tests generated" << endl;
        ENDL(3);
    }

    //********************************************************************************************

    if(generate_nonuniform_graphs){
        int test_id = 1;
        clog << "Generating nonuniform tests" << endl;

        VI nodes = {
                2'000, 5'000, 10'000, 20'000, 50'000
        };

        VI degrees = {
                10,20,50,100
        };


        VPII random_sizes = StandardUtils::productIf(nodes, degrees, [&](int a, int b){
            return a*b <= MAX_E;
        });

        for( auto & [N,M] : random_sizes ) M *= N;

        for( auto [N,M] : random_sizes ){
            VI probabs(N);

            {
                iota(ALL(probabs), N);
                clog << "\rCreating test #" << test_id << flush;
                ofstream str(out_directory + sep + "F_" + getTestName(test_id) + ".txt");
                createRandomNonuniform(str, N, M, probabs);
                str.close();
                test_id++;
            }

            {
                iota(ALL(probabs), sqrt(N));
                clog << "\rCreating test #" << test_id << flush;
                ofstream str(out_directory + sep + "F_" + getTestName(test_id) + ".txt");
                createRandomNonuniform(str, N, M, probabs);
                str.close();
                test_id++;
            }
        }

        clog << endl << "nonunifrom tests generated" << endl;
        ENDL(3);
    }


    clog << "Tests generated!" << endl;
}


void CETestGraphGenerator::writeGraphToStream(ostream &str, VVI V, string message) {
    assert(GraphUtils::isSimple(V));
    VPII edges = GraphUtils::getGraphEdges(V);
    int N = V.size();
    int M = edges.size();

    str << "c " << message << endl;
    str << "p ce " << N << " " << M << endl;
    for(auto [a,b] : edges) str << a+1 << " " << b+1 << endl;
}


void CETestGraphGenerator::createRandom1(ostream &str, int N, double p) {
    VVI V = GraphGenerator::getRandomGraph(N,p);
    assert(GraphUtils::isSimple(V));

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "random graph N: " << N << ", p: " << p;

    writeGraphToStream(str,V, s.str());
}

void CETestGraphGenerator::createRandom2(ostream &str, int N, int M) {
    VVI V = GraphGenerator::getRandomGraph(N,M);
    assert(GraphUtils::isSimple(V));
    assert( M == GraphUtils::countEdges(V) );

    double dns = 2.0*M / ( 1ll * N * (N-1) );

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "random graph N: " << N << ", M: " << M << ", dns: " << dns;
    writeGraphToStream(str,V, s.str());
}

void CETestGraphGenerator::createRandomKPartite(ostream &str, int M, int K, double p, int minC, int maxC) {
    UniformIntGenerator rnd(minC, maxC);
    VI partition_sizes(K);
    for( int i=0; i<K; i++ ) partition_sizes[i] = rnd.rand();
    VVI V = GraphGenerator::getRandomKPartiteGraph( partition_sizes, M );
    assert(GraphUtils::isSimple(V));
    int N = V.size();
//    assert( M == GraphUtils::countEdges(V) );

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "random " << K << "-partite graph N: " << N << ", M: " << M << ", p: " << p;
    s << ", minC: " << minC << ", maxC: " << maxC;
//    s << ", partition sizes: "; for(int c : partition_sizes) s << c << " ";

    writeGraphToStream(str,V, s.str());
}

void CETestGraphGenerator::createRandomKPartite(ostream &str, int M, double p, VI partition_sizes) {
    VVI V = GraphGenerator::getRandomKPartiteGraph( partition_sizes, M );
    assert(GraphUtils::isSimple(V));
    int N = V.size();
    int K = partition_sizes.size();
    double dns = 2.0*M / ( 1ll * N * (N-1) );

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "random " << K << "-partite graph N: " << N << ", M: " << M << ", p: " << p << ", dns: " << dns;
    if( partition_sizes.size() <= 10 ){
        s << ", partition sizes: "; for(int c : partition_sizes) s << c << " ";
    }
    else s << ", minC: " << *min_element(ALL(partition_sizes)) << ", maxC: " << *max_element(ALL(partition_sizes));

    writeGraphToStream(str,V, s.str());

}

void CETestGraphGenerator::createRandomRegular(ostream &str, int N, int K) {
    assert( (N*K)%2 == 0 );

    int iters = min( 1ll*N*N/2, 5'000'000ll );
    VVI V = GraphGenerator::getRandomKRegularGraph(N,K, iters);
    assert(GraphUtils::isSimple(V));
    int M = (N*K/2);
    assert( M == GraphUtils::countEdges(V) );
    double dns = 2.0*M / ( 1ll * N * (N-1) );

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "random " << K << "-regular graph N: " << N << ", M: " << M << ", dns: " << dns;
    writeGraphToStream(str, V, s.str());
}

void CETestGraphGenerator::createRandomClustered1(ostream &str, int K, double p, double q, int minC, int maxC) {
    UniformIntGenerator rnd(minC, maxC);
    VI cluster_sizes(K);
    for( int i=0; i<K; i++ ) cluster_sizes[i] = rnd.rand();

    VVI V = GraphGenerator::getRandomClusteredGraph1( cluster_sizes, p, q );
    assert(GraphUtils::isSimple(V));
    int N = V.size();
    int M = GraphUtils::countEdges(V);

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "clustered graph with " << K << " clusters, N: " << N << ", M: " << M;
    s << ", p: " << p << ", q: " << q << ", minC: " << minC << ", maxC: " << maxC;
//    s << ", cluster sizes: "; for(int c : cluster_sizes) s << c << " ";

    writeGraphToStream(str,V, s.str());
}

void CETestGraphGenerator::createRandomClustered1(ostream &str, VI cl_sizes, double p, double q, int UB, int E) {
    VVI V = GraphGenerator::getRandomClusteredGraph1( cl_sizes, p, q );

    if(E != -1){
        int realE = GraphUtils::countEdges(V);
//        DEBUG(E);
//        DEBUG(realE);
        assert( E == realE );
    }

    assert(GraphUtils::isSimple(V));
    int N = V.size();
    int M = GraphUtils::countEdges(V);
    int K = cl_sizes.size();

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "clustered graph with " << K << " clusters, N: " << N << ", M: " << M;
    s << ", p: " << p << ", q: " << q << ", minC: " << cl_sizes[0] << ", maxC: " << cl_sizes.back();
    if(UB != -1) s << ", UB: " << UB;
//    s << ", cluster sizes: "; for(int c : cluster_sizes) s << c << " ";

    writeGraphToStream(str,V, s.str());
}


void CETestGraphGenerator::createRandomNonuniform(ostream &str, int N, int M, VI probabs) {
    VVI V = GraphGenerator::getRandomGraphNonuniform(N, M, probabs);
    assert(GraphUtils::isSimple(V));
    double dns = 2.0*M / ( 1ll * N * (N-1) );

    stringstream s;
    s << fixed;
    s.precision(4);
//    s << "random graph with nonuniform edge selection probabilities, N: " << N << ", M: " << M;
    s << "random, nonuniform, N: " << N << ", M: " << M << ", dns: " << dns;
    {
        int min_pr = *min_element(ALL(probabs));
        int max_pr = *max_element(ALL(probabs));
        s << ", min_pr: " << min_pr << ", max_pr: " << max_pr;
    }
    writeGraphToStream(str,V, s.str());
}

void CETestGraphGenerator::createRandomSpatialEuclid(ostream & str, int N, int k, int rand_neighs, double sparsity) {
    int U = 1e7;
    using namespace StandardUtils;
    using namespace CombinatoricUtils;
    VPII points = zip( getRandomSubset(U,N), getRandomSubset(U,N) );
    VVI V = GraphGenerator::getSpatialGraph( N, k, points, GraphGenerator::getEuclideanMetric(), rand_neighs, sparsity );
    assert(GraphUtils::isSimple(V));
    int M = GraphUtils::countEdges(V);

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "spatial, euclidean metric, N: " << N << ", M: " << M << ", k: " << k
//      << ", rand_neighs: " << rand_neighs
      << ", sparsity: " << sparsity;
    writeGraphToStream(str,V, s.str());
}

void CETestGraphGenerator::createRandomSpatialEuclidTorus(ostream & str, int N, int k, int rand_neighs, double sparsity) {
    int U = 1e7;
    using namespace StandardUtils;
    using namespace CombinatoricUtils;
    VPII points = zip( getRandomSubset(U,N), getRandomSubset(U,N) );
    VVI V = GraphGenerator::getSpatialGraph( N, k, points,GraphGenerator::getEuclideanMetricOnTorus(0,U,0,U), rand_neighs, sparsity );
    assert(GraphUtils::isSimple(V));
    int M = GraphUtils::countEdges(V);

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "spatial, euclidean metric on torus, N: " << N << ", M: " << M << ", k: " << k
//      << ", rand_neighs: " << rand_neighs
      << ", sparsity: " << sparsity;
    writeGraphToStream(str,V, s.str());
}

void CETestGraphGenerator::createRandomSpatialManhattan(ostream & str, int N, int k, int rand_neighs, double sparsity) {
    int U = 1e7;
    using namespace StandardUtils;
    using namespace CombinatoricUtils;
    VPII points = zip( getRandomSubset(U,N), getRandomSubset(U,N) );
    VVI V = GraphGenerator::getSpatialGraph( N, k, points, GraphGenerator::getManhattanMetric(), rand_neighs, sparsity );
    assert(GraphUtils::isSimple(V));
    int M = GraphUtils::countEdges(V);

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "random spatial graph with manhattan metric, N: " << N << ", M: " << M << ", k: " << k
      << ", rand_neighs: " << rand_neighs;
    writeGraphToStream(str,V, s.str());
}

void CETestGraphGenerator::createRandomSpatialManhattanTorus(ostream & str, int N, int k, int rand_neighs, double sparsity) {
    int U = 1e7;
    using namespace StandardUtils;
    using namespace CombinatoricUtils;
    VPII points = zip( getRandomSubset(U,N), getRandomSubset(U,N) );
    VVI V = GraphGenerator::getSpatialGraph( N, k, points, GraphGenerator::getManhattanMetricOnTours(0,U,0,U), rand_neighs, sparsity );
    assert(GraphUtils::isSimple(V));
    int M = GraphUtils::countEdges(V);

    stringstream s;
    s << fixed;
    s.precision(4);
    s << "random spatial graph with manhattan metric on torus, N: " << N << ", M: " << M << ", k: " << k
      << ", rand_neighs: " << rand_neighs;
    writeGraphToStream(str,V, s.str());
}




