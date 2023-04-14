//
// Created by sylwester on 9/6/19.
//

#ifndef ALGORITHMSPROJECT_STANDARDUTILS_H
#define ALGORITHMSPROJECT_STANDARDUTILS_H

#include "Makros.h"
#include "RandomNumberGenerators.h"

namespace StandardUtils{

    /**
     * Equivalent to python slicing. Extracts and returns in a vector all elements from [a,b) that are at positions
     * a, a+step, a+2*step, ...
     * It is also possible to pass a > b and step < 0
     */
    template<class _T>
    vector<_T> slice( vector<_T> & v, int a, int b, int step=1 ){
        assert(step != 0);
        vector<_T> res;
        res.reserve( abs( (b-a) / step ) );
        if(a<=b){
            assert(step > 0);
            while(a<b){
                res.push_back( v[a] );
                a += step;
            }
        }
        else{
            assert(step < 0);
            while(a>b){
                res.push_back(v[a]);
                a += step;
            }
        }
        return res;
    }



    template<class _T, class _rnd>
    void shuffle( vector<_T> & V, _rnd rnd ){
        std::uniform_int_distribution<long long> unif( 0, 10ll * V.size() );

        for( int i=(int)V.size()-1; i>=0; i-- ){
            int ind = unif(rnd) % (i+1);
            if( ind != i ) swap( V[i], V[ind] );
        }
    }

    template<class _T>
    void shuffle( vector<_T> & V ){
        UniformIntGenerator rnd(0,1e9);
        shuffle(V, rnd.getRNG());
    }

    /**
     * Creates a vector of pairs from two vectors. E.g. from [0,4,7], ['a', 's'] it will make [ {0,'a'}, {4,'s'} ]
     */
    template<class _t, class _s>
    vector< pair<_t,_s> > zip( vector<_t> a, vector<_s> b ){
        vector<pair<_t,_s>> res(min(a.size(), b.size()));
        for( int i=0; i<min(a.size(), b.size()); i++ ) res[i] = {a[i], b[i]};
        return res;
    }

    /**
     * Creates all possible pairs (x,y) with x \in a and y \in b.
     */
    template<class T>
    vector<pair<T,T>> product( vector<T> a, vector<T> b ){
        vector<pair<T,T>> res;
        res.reserve(a.size() * b.size());

        for( int i=0; i<a.size(); i++ ){
            for(int j=0; j<b.size(); j++){
                res.emplace_back( a[i], b[j] );
            }
        }
        return res;
    }

    /**
     * Does the product if the predicate [pred] returns true.
     * [pred] must be callable pred(_T,_T) and return a boolean value
     */
    template<class T, class _P>
    vector<pair<T,T>> productIf( vector<T> a, vector<T> b, _P pred ){
        vector<pair<T,T>> res;
        res.reserve(a.size() * b.size());

        for( int i=0; i<a.size(); i++ ){
            for(int j=0; j<b.size(); j++){
                if( pred(a[i], b[j]) ) res.emplace_back( a[i], b[j] );
            }
        }
        return res;
    }

}

#endif //ALGORITHMSPROJECT_STANDARDUTILS_H
