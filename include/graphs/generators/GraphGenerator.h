//
// Created by sylwester on 8/8/19.
//

#ifndef ALGORITHMSPROJECT_GRAPHGENERATOR_H
#define ALGORITHMSPROJECT_GRAPHGENERATOR_H

#include "Makros.h"

class GraphGenerator {

public:

    /**
     * Creates graph G(N,M)
     * @param N
     * @param M
     * @param directed if true, then a directed graph is created, otherwise undirected
     * @return random graph with N vertices and M edges.
     */
    static VVI getRandomGraph(int N, int M, bool directed = false);

    /**
     * Creates random graph G(N,p)
     * @param N
     * @param p
     * @return random graph with N vertices, in which every edge hash probability p of being in that graph
     */
    static VVI getRandomGraph(int N, double p, bool directed = false);

    /**
     *
     * @param L
     * @param R
     * @param M
     * @return random bipartite graph with M edges, where the first L vertices are in first part, and last R are in second part of the bipartition
     */
    static VVI getRandomBipartiteGraph( int L, int R, int M );

    /**
     * Creates a random K-partite graph with M randomly chosen edges, where K = partition_sizes.size(). i-th partition of the graph will have
     * partition_sizes[i] nodes (mutually disjoint).
     */
    static VVI getRandomKPartiteGraph( VI partition_sizes, int M );

    /**
     * Creates a random graph with N nodes and M edges. Each edges is chosen as a pair of nodes. Nodes are selected
     * randomly with probability of choosing node v equal to probabs[v] / S, where S is the sum of elements in probabs.
     *
     * #CAUTION!! Edges are selected each time randomly. Expected time to create the graph is therefore
     * at least M * logM (according to the Coupon-collector problem it will be M*logM if probabilities are equal),
     * but it may MUCH MUCH SLOWER, if the probabilities have a very high variance and the graph is dense.
     */
    static VVI getRandomGraphNonuniform(int N, int M, VI probabs, bool directed = false);

    /**
     *
     * @param L
     * @param R
     * @param p
     * @param directed
     * @return random bipartite graph in which each edge is with probability p, where the first L vertices are in first part, and last R are in second part of the bipartition
     */
    static VVI getRandomBipartiteGraph( int L, int R, double p );

    /**
     * Creates and returns random grid with rows*columns nodes. Node in row r and column c will have id = r*rows + c
     * @param rows
     * @param columns
     * @return random grid with @{rows} rows and @{columns} columns.
     */
    static VVI getGrid(int rows, int columns );

    /**
     * Function generates random tree. In first step it generates random prufer code, then creates tree represented by that code.
     * @param N
     * @return random tree.
     */
    static VVI getRandomTreePrufer(int N);

    /**
     * Function creates and returns path on N vertices
     * @param N
     * @return path on N vertices
     */
    static VVI getPath(int N, bool randomOrder = true);

    /**
     * Creates a graph by selecting randomly [triangles] triangles (a,b,c), then adding all edges of
     * those triangles to the graph
     *
     * @param directed if true, then directed edges will be added, otherwise undirected edges
     */
    static VVI getUnionOfTriangles(int N, int triangles, bool directed = false);

    /**
     * @return a start with N nodes, node 0 is the center of the star
     */
    static VVI getStar( int N );

    /**
     * Craetes and returns a full graph K_n
     */
    static VVI getClique(int N);

    /**
     * Creates a random k-regular graph.
     * #CAUTION! The condition for existance of k-regular graphs needs to be met, t
     * hat is N >= K+1 and  N*K == 0 (mod 2)
     *
     * If K is even, then we create an initial graph as a 'circulant' graph, by adding to each node i neighbors
     * (i+j) % N, for j =1,2,...,K/2.
     *
     * If K is odd, then N must be even. If K is small, then we create randomly K full matchings, then take a union
     * of those matchings. A matching is created by taking a random permutation of nodes, then taking edges as two
     * consecutive nodes.
     * If, however, K is large, then we first craete a clique K_N, then K times find a perfect matching in that graph,
     * and remove that matching. This way we will always find a K edge-disjoint matchings.
     *
     * After creating an initial graph, we consider [iters] (or N*N if [iters] is -1) pairs of present edges,
     * (a,b) and (c,d), then we change that pair to (a,c) and (b,d).
     *
     * @param N number of nodes
     * @param K degree of each vertex
     */
    static VVI getRandomKRegularGraph( int N, int K, int iters = -1 );

    /**
     * Expands the graph in the following way:
     * For each node v in V, two copies v1 and v2 will be creates.
     * Then for each edge (a,b) in V, edges (a1,b1), (a2,b2), (a1,b2), (a2,b1) will be in the resulting graph.
     */
    static VVI expand1( VVI V );

    /**
     * Expands the graph in the following way:
     * For each node v in V, two copies v1 and v2 will be creates.
     * Then for each edge (a,b) in V, edges (a1,b2), (a2,b1) will be in the resulting graph.
     */
    static VVI expand2( VVI V );

    /**
     * Expands the graph in the following way:
     * For each node v in V, two copies v1 and v2 will be creates.
     * Then for each edge (a,b) in V, edges (a1,b1), (a2,b2) will be in the resulting graph.
     * Additionally (v1,v2) will be in the graph for each v.
     */
    static VVI expand3( VVI V );

    /**
     * Expands the graph in the following way:
     * For each node v in V, two copies v1 and v2 will be creates.
     * Then for each edge (a,b) in V, edges (a1,b1), (a2,b2), (a1,b2), (a2,b1) will be in the resulting graph.
     * Additionally (v1,v2) will be in the graph for each v.
     */
    static VVI expand4( VVI V );



    /**
     * Creates an isomorphic graph, by taking a random permutation, then changing ids.
     */
    static VVI createIsomorphic(VVI V);

    /**
     * Makes a union of two given graphs.
     * In the first step the smaller graph with n nodes will be mapped to the larger graph with N >= n nodes. This is
     * done by remapping nodes using a random injection [n] -> [N].
     *
     * Then, in the resulting graph, an edge will be present if it is present in either of two graphs (the larger one,
     * or the smaller one mapped to the larger size).
     */
    static VVI mergeGraphs1( VVI V, VVI H );

    /**
     * Creates a graph by taking N points on a plane with integer coordinates.
     * -10^9 \leq \leq x,y \leq 10^9. Then, for each node v, chooses k nodes that are closest to this node and edge
     * between node v and those closest nodes will be added to the graph.
     * Note, that edges may repeat for different v.
     *
     * Uses function [dist] to calculate the distance between two points.
     *
     * Generating takes time at least O( N^2 * log(N) ) - for each node v we sort all other nodes using [dist] function
     *
     * If [rand_neighs] > 0, the from closest k neighbors of each node we will select random [rand_neighs] nodes to
     * add connections between. This may be used to try to 'sparsify' clique that usually are created if we add all
     * connections for closest neighborhoods
     *
     * After the graph is created, only [sparsity] percentage of randomly chosen edges will be retained in the graph.
     */
    static VVI getSpatialGraph(int N, int k, VPII coords, function<double(PII, PII)> dist, int rand_neighs = 0,
                               double sparsity = 1.0 );

    /**
     * The following functions create and return functions calculating distances between two points (x1,y1), (x2,y2).
     * In case of calculating distances on Torus, the bounding rectangle needs to be provided as well.
     */
    static function<double(PII, PII)> getEuclideanMetric();
    static function<double(PII, PII)> getEuclideanMetricOnTorus(int minx, int maxx, int miny, int maxy);
    static function<double(PII, PII)> getManhattanMetric();
    static function<double(PII, PII)> getManhattanMetricOnTours(int minx, int maxx, int miny, int maxy);

    /**
     * Transforms given [metric] to work on torus. This is done by considering for a pair (a,b) 5 points:
     * (a,b), (a+W,b), (a-W,b), (a,b+H), (a,b-H)
     * @param metric metric for which corresponding metric on torus should be returned
     */
    static function<double(PII, PII)> getMetricForTorus(int minx, int maxx, int miny, int maxy, function<double(PII, PII)> metric);


    /**
     * Creates a random 'clustered' graph.
     * That is, there will be cluster_sizes.size() 'clusters'.
     * An edge (a,b) with nodes in the same cluster will 'have probability p' - for a cluster of size C, exactly
     * p*C*(C-1)/2 random edges will be selected.
     *
     * An edge with ends in different clusters will 'have probability q' - if there are X such possible edges, then
     * q*X of them will be chosen randomly.
     *
     * By selecting edges this way, we ensure that there are no 'large' deviations for small graphs/clusters and
     * that time complexity remains roughly O(E), where E is the total number of edges, hence large sparse graphs can
     * be created this way.
     *
     * #CAUTION!! Current implementation works in time O(N^2) !!
     */
    static VVI getRandomClusteredGraph1(VI cluster_sizes, double p, double q );

};


#endif //ALGORITHMSPROJECT_GRAPHGENERATOR_H
