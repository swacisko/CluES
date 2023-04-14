//
// Created by sylwester on 1/17/22.
//

#ifndef ALGORITHMSPROJECT_CETESTGRAPHGENERATOR_H
#define ALGORITHMSPROJECT_CETESTGRAPHGENERATOR_H

#include "Makros.h"

class CETestGraphGenerator{

public:

    string out_directory = "ClusterEditingTests";

    /**
     * Generates a collection of instances for CE.
     */
    void generateAllInstances();

    /**
     * Other version of generating instances - based on densities instead of degrees, and some merged categories
     */
    void generateAllInstances2();

    void generateAllInstances3();

    /**
     * Final test for post-PACE publication. Generate 10 instances of each graph class.
     */
    void generateAllInstances4();

    /**
     * Writes given graph in DIMACS format to given stream.
     * @param str
     * @param V
     */
    void writeGraphToStream(ostream & str, VVI V, string message);

    /**
     * Creates and write tha random graph.
     * Graph are generated with N vertices and random p * N*(N-1)/2 edges, where p is a probability of a single edge.
     * Erdos-Renyi model.
     */
    void createRandom1(ostream& str, int N, double p);


    /**
     * Creates and write tha random graph.
     * Graph are generated with N vertices and random M edges to the stream [str].
     */
    void createRandom2(ostream& str, int N, int M);

    /**
     * Creates a random K-partite graph with M edges. Each of K partition sets cluster has a random size from
     * set [minC, maxC]. The number of nodes is therefore not know from start, unless minC == maxC.
     * Value of [p] is used only in comment line about the graph generated.
     */
    void createRandomKPartite(ostream & str,  int M, int K, double p, int minC, int maxC);
    void createRandomKPartite(ostream & str,  int M, double p, VI partition_sizes );

    /**
     * Creates and write tha random K-regular graph with N vertices and M edges to the stream [str].
     * #CAUTION!
     * Value N*K must be even and N > K
     */
    void createRandomRegular(ostream& str, int N, int K);


    /**
     * Creates a random 'clustered' graph.
     * That is, there will be K 'clusters'.
     * Graph will be generated using GraphGenerator::getRandomClusteredGraph1
     *
     * Clusters will have random sizes from set [minC, maxC], so the number of nodes N is not know in advance.
     */
    void createRandomClustered1( ostream & str, int K, double p, double q, int minC, int maxC );

    /**
     *
     * @param str
     * @param cl_sizes
     * @param p probability of edge inside clusters
     * @param q probability of edges between differenct clusters
     * @param UB upper bound on solution size (by making given clusters clusters :) )
     * @param E expected number of edges - for assertion of UB correctness
     */
    void createRandomClustered1(ostream &str, VI cl_sizes, double p, double q, int UB, int E);

    /**
     * Uses GraphGenerator::getRandomGraphNonuniform() to generate the graph
     */
    void createRandomNonuniform( ostream & str, int N, int M,  VI probabs );

    /**
     * Creates a random tree from a random Prufer code, using GraphGenerator::getRandomTreePrufer() function,
     * then adds [additional_edges] random edges that are not present in the tree.
     */
    void createRandomTree( ostream & str, int N, int additional_edges);

    /**
     * Creates a graph byu generating [trees] random trees, then taking a graph that contains all edges in those trees.
     */
    void createRandomTreeUnion(  ostream & str, int N, int trees );


    /**
     * Creates a grid with N rows and M edges, then adds [additional_edges] additional edges that were not present in
     * the mesh.
     */
    void createRandomGrid( ostream & str, int N, int M, int additional_edges);


    /**
     * The same as [createRandomSpatial1], but we take Euclidean distance.
     */
    void createRandomSpatialEuclid( ostream & str, int N, int k, int rand_neighs, double sparsity );

    /**
    * The same as [createRandomSpatial1], but we take Euclidean distance on a Torus.
    */
    void createRandomSpatialEuclidTorus(ostream & str, int N, int k, int rand_neighs, double sparsity );

    /**
     * The same as [createRandomSpatial1], but we take Manhattan distance.
     */
    void createRandomSpatialManhattan(ostream & str, int N, int k, int rand_neighs, double sparsity );

    /**
    * The same as [createRandomSpatial1], but we take Manhattan distance on a Torus.
    */
    void createRandomSpatialManhattanTorus(ostream & str, int N, int k, int rand_neighs, double sparsity );

};

#endif //ALGORITHMSPROJECT_CETESTGRAPHGENERATOR_H
