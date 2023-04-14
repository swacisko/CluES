//
// Created by sylwester on 8/8/19.
//

#ifndef ALGORITHMSPROJECT_GRAPHINDUCER_H
#define ALGORITHMSPROJECT_GRAPHINDUCER_H

#include "Makros.h"


struct InducedGraph{
    VVI *par; // parent graph;
    VI nodes; // this is a vector of vertices that induce the graph.
    // vertex with number nodes[i] has in induced graph number i.

    unordered_map<int,int> perm; // perm[t] is the number d such that nodes[d] = t; E.g. if graph is induced by [2,8,5] then perm[8] = 1.  So perm[ nodes[i] ] = i for i in [ 0,SIZE(V) ) and
    // nodes[ perm[i] ] = i for i in {nodes[0], nodes[1], ..., nodes.back() }

    VPII edges; // this is a vector of edges that induce a graph. It may be empty if the graph is induced by nodes
    VVI V; // induced graph

    friend ostream& operator<<(ostream& str, InducedGraph& g);

    void write(){
        cerr << "Graph induced by: " << flush; WRITE(nodes);
        if( !edges.empty() ){
            cerr << "induced by edges:" << endl;
            REP(i,SIZE(edges)) cerr << WRP(edges[i]) << endl;
        }
        WRITE_ALL( V, "Graph structure",0 );

    }
};

/**
 * For weighted graphs on structure VVPII
 */
struct InducedGraphPI{
    VVPII *par; // parent graph;
    VI nodes; // this is a vector of vertices that induce the graph.
    // vertex with number nodes[i] has in induced graph number i.

    unordered_map<int,int> perm; // perm[t] is the number d such that nodes[d] = t; E.g. if graph is induced by [2,8,5] then perm[8] = 1.  So perm[ nodes[i] ] = i for i in [ 0,SIZE(V) ) and
    // nodes[ perm[i] ] = i for i in {nodes[0], nodes[1], ..., nodes.back() }

    VPII edges; // this is a vector of edges that induce a graph. It may be empty if the graph is induced by nodes
    VVPII V; // induced graph

    friend ostream& operator<<(ostream& str, InducedGraphPI& g);
};


class GraphInducer{
public:

    // return graph induced by given nodes. Works for directed graphs as well (V can be directed).
    static InducedGraph induce( VVI & V, VI & nodes );

    /**
     * For weighted graphs on structure VVPII
     */
    static InducedGraphPI induce(VVPII & V, VI & nodes );

    // returns graph induced by given edges. Works for directed graphs (V can be directed) as welll
    // if directed == true then each edge in edges will be treated as directed edge. Otherwise it will be treated as undirected, bidirectional edge.
    static InducedGraph induce( VVI & V, VPII & edges, bool directed = false );

};

#endif //ALGORITHMSPROJECT_GRAPHINDUCER_H
