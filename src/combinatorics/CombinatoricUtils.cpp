//
// Created by sylwester on 8/27/19.
//

#include "combinatorics/CombinatoricUtils.h"
#include <utils/RandomNumberGenerators.h>
#include <combinatorics/CombinatoricUtils.h>


namespace CombinatoricUtils{


    VI getRandomPermutation(int N){
        UniformIntGenerator rnd(0,1e9);
        return getRandomPermutation(N, rnd.rand());
//        VI perm(N);
//        iota(ALL(perm),0);
//        random_shuffle(ALL(perm));
//        return perm;
    }

    VI getRandomPermutation(int N, unsigned seed){
        VI perm(N);
        iota(ALL(perm),0);
        shuffle(ALL(perm), mt19937(seed) );
        return perm;
    }


    VI getRandomSequence( int U, int N ) {
        VI seq(N);
        UniformIntGenerator g( 0,U );
        for(int i=0; i<N; i++) seq[i] = g.rand();
        return seq;
    }


    VI getRandomSubset( int U, int L ){
        UniformIntGenerator rnd(0,1e9);
        return getRandomSubset(U,L,rnd.rand());
    }

    VI getRandomSubset( int U, int L, unsigned seed ){
        UniformIntGenerator gen(0,U, seed);

        VI res;
        if( L < U / 20 ){
            unordered_set<int> zb;
            while( zb.size() < L ) zb.insert( gen.rand() );
            res = VI(ALL(zb));
        }else{
            VI perm = getRandomPermutation(U+1, seed);
            if( L < perm.size() ) return VI( perm.begin(), perm.begin() + L );
            else return VI( ALL(perm) );
        }
        return res;
    }

    VLL getRandomSubset( LL U, LL L, unsigned seed ){
        assert( L <= U );
        UniformIntGenerator gen(0,U, seed == 1 ? RandomNumberGenerators::DEFAULT_SEED : seed);

        VLL res;
        if( L < U / 20 ){
            unordered_set<LL> zb;
            while( zb.size() < L ) zb.insert( gen.rand() );
            res = VLL(ALL(zb));
        }else{
            VI permint = getRandomPermutation(U+1, seed);
            res.resize(L);
            copy( permint.begin(), permint.begin()+L, res.begin() );
        }
        return res;
    }


    VI getFullSetDifference(int N, VI A) {
        VB helper(N,false);
        for(int d : A) helper[d] = true;
        VI res; res.reserve(N-A.size());
        for(int i=0; i<N; i++) if(!helper[i]) res.push_back(i);
        return res;
    }

    void allPartitions(VI &sets, function<void(VI &, int ind, int d, bool added)> fun) {

        VI seq;
        function< void( int ind ) > generator = [&generator, &sets, &fun, &seq]( int ind ){

                for( int i=0; i<sets[ind]; i++ ){
                    seq.push_back(i);
                    fun( seq, ind,i,true );

                    if( ind +1 < sets.size() ) generator(ind+1);

                    seq.pop_back();
                    fun( seq, ind,i,false );

                }
        };

        generator(0);
    }


}

