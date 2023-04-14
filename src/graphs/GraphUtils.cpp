//
// Created by sylwester on 8/8/19.
//


#include <graphs/GraphUtils.h>

VI GraphUtils::getComplimentaryNodes( VVI & V, VI & nodes ){
    VB inNodes( V.size(),false );
    for(auto p : nodes) inNodes[p] = true;
    VI res;
    for( int i=0; i<V.size(); i++ ){
        if( !inNodes[i] ) res.push_back(i);
    }
    return res;
}


VPII GraphUtils::getGraphEdges( VVI & V, bool directed ){
    if(directed) return getDirectedGraphEdges(V);

    VPII res;
    res.reserve( countEdges(V) );
    for( int i=0; i<V.size(); i++ ){
        for( int d : V[i] ){
            if(d>i) res.push_back( {i,d} );
        }
    }
    return res;
}

VPII GraphUtils::getDirectedGraphEdges( VVI & V ){
    VPII res;
    res.reserve( countEdges(V) );
    for( int i=0; i<V.size(); i++ ){
        for( int d : V[i] ){
            res.push_back( {i,d} );
        }
    }
    return res;
}

VI GraphUtils::getNeighborhood(VVI &V, VI &S, bool useSet) {
    VI res;
    if( useSet ){
        unordered_set<int> zb;
        for( int p : S ){
            for(int d : V[p]){
                zb.insert(d);
            }
        }
        res = VI(ALL(zb));
    }else{
        VB was(V.size(),false);
        for( int p : S ){
            for(int d : V[p]){
                if( !was[d] ){
                    was[d] = true;
                    res.push_back(d);
                }
            }
        }
    }
    return res;
}


VVI GraphUtils::transposeGraph(VVI &v) {
    VVI g( SIZE(v) );
    REP( i,SIZE(v) )	REP( k,SIZE(v[i]) )	g[ v[i][k] ].PB(i);
    return g;
}


void GraphUtils::addEdge(VVI &V, int a, int b, bool directed) {
    V[a].push_back(b);
    if( !directed ) V[b].push_back(a);
}

void GraphUtils::removeEdge(VVI &V, int a, int b, bool directed) {
    auto rem = [ &V ](int a, int b){
        for( int i=(int)V[a].size()-1; i>=0; i-- ){
            if( V[a][i] == b ){
                swap(V[a][i], V[a].back() );
                V[a].pop_back();
                break;
            }
        }
    };


    rem(a,b);
    if(!directed) rem(b,a);

}


int GraphUtils::countEdges(VVI &V, const bool directed) {
    int res = 0;
    for(auto& v : V) res += v.size();
    if(!directed) return res >> 1;
    else return res;
}

VVI GraphUtils::sortNodeNeighborhoods( VVI & V ){
    VVI V2(V.size());
    for(int i=0; i<V.size(); i++) V2[i].reserve(V[i].size());
    for( int i=0; i<V.size(); i++ ){
        for( int d : V[i] ) V2[d].push_back(i);
    }
    return V2;
}

void GraphUtils::removeNodesFromGraph(VVI &V, VI nodes) {

    unordered_set<int> nd(ALL(nodes));
    VI neigh = getNeighborhood( V,nodes,true );
    for( int t : neigh ){
        if( nd.count(t) ) continue;

        for( int i = (int)V[t].size()-1; i>=0; i-- ){
            if( nd.count(V[t][i]) ){
                swap( V[t][i], V[t].back() );
                V[t].pop_back();
            }
        }
    }

    for( int d : nodes ) V[d].clear();
}

void GraphUtils::removeNodesFromGraph(VVI &V, VVI &W, VI nodes) {
    unordered_set<int> nd(ALL(nodes));
    VI neigh = getNeighborhood( V,nodes,true );
    for( int t : neigh ){
        if( nd.count(t) ) continue;

        for( int i = (int)V[t].size()-1; i>=0; i-- ){
            if( nd.count(V[t][i]) ){
                swap( V[t][i], V[t].back() );
                V[t].pop_back();

                swap( W[t][i], W[t].back() );
                W[t].pop_back();
            }
        }
    }

    for( int d : nodes ){
        V[d].clear();
        W[d].clear();
    }
}

void GraphUtils::writeGraphHumanReadable(VVI &V) {
    clog << "****" << endl;
    for( int i=0; i<V.size(); i++ ){
        clog << i << ": ";
        VI neigh = V[i];
        sort(ALL(neigh));
        for(int d : neigh) clog << d << "  ";
        clog << endl;
    }
}


void GraphUtils::removeEdges(VVI &V, VPII &edges, bool directed) {

    if( directed == false ){
        int E = edges.size();
        for( int i=0; i<E; i++ ) edges.emplace_back( edges[i].second, edges[i].first ); // adding reverse edges to remove
    }

    sort( ALL(edges) );

    for( int i=0; i<edges.size(); i++ ){
        int p = i;
        unordered_set<int> toRemove;
        while( p < edges.size() && edges[p].first == edges[i].first ){
            toRemove.insert( edges[p].second );
            p++;
        }

        int t = edges[i].first;
        for( int k=(int)V[t].size()-1; k>=0; k-- ){
            int d = V[t][k];
            if( toRemove.count(d) ){
                swap( V[t][k], V[t].back() );
                V[t].pop_back();
            }
        }

        i = p-1;
    }

}

bool GraphUtils::isConnected(VVI &V) {
    VB was(V.size(),false);
    int cnt = 0;
    function< void(int) > dfs = [&V,&was, &dfs, &cnt](int num){
        was[num] = true;
        cnt++;
        for( int d : V[num] ) if( !was[d] ) dfs(d);
    };
    dfs(0);
    return (cnt == V.size());
}


int GraphUtils::regularity(VVI V) {
    int k = V[0].size();
    for(int i=1; i<V.size(); i++) if( V[i].size() != k ) return -1;
    return k;
}

VVI GraphUtils::getGraphForEdges(VPII edges, bool directed) {
    int N = 0;
    for(auto & [a,b] : edges) N = max(N, max(a,b));
    VVI V(N+1);
    for( auto & [a,b] : edges ) addEdge(V,a,b,directed);
    return V;
}

VVI GraphUtils::remapGraph(VVI V, VI perm) {
    int N = V.size();
    VVI H(N);
    for( int i=0; i<N; i++ ){
        H[perm[i]].reserve(V[i].size());
        for(int d : V[i]) H[perm[i]].push_back(perm[d]);
    }

    return H;
}

bool GraphUtils::isSimple(VVI V) {
    int N = V.size();
     VB helper(N,false);

     for( int i=0; i<N; i++ ){
         for( int d : V[i] ){
             if( d == i || helper[d] ) return false;
             helper[d] = true;
         }

         for(int d : V[i]) helper[d] = false;
     }

     return true;
}


