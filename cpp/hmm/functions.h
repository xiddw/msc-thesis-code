#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>

#include <utility>

typedef unsigned int uint;

// Normalize a boost::vector to sum up to one
template<typename V> 
void normalize_row(V &v) { 
   v /= sum(v); 
} 

// Normalize a boost::matrix in order to have its rows/columns to sum up to one
template<typename M> 
void normalize(M &m, bool byrow = true) { 
  using namespace boost::numeric::ublas;

  for(uint i = 0; i<m.size1(); ++i) { 
     auto r = row(m, i);
     normalize_row(r);
   }
} 

// Resize two vectors to have the same lenght repeating and/or truncating some values.
template<typename D>
std::pair< vector<D>, vector<D> > resize(vector<D> a, vector<D> b) { 
  bool alessthanb = a.size() < b.size();

  vector<D> aa = alessthanb ? a : b; 
  vector<D> bb = alessthanb ? b : a;

  uint la = aa.size();
  uint lb = bb.size();

  uint nr = uint(lb / la);

  vector<D> tt = zero_vector<D>(nr *la);

  for(int i=0; i<la; ++i) {
    for(int r=0; r<nr; ++r) {
      tt(i*nr + r) = aa(i);
    }
  }

  uint lt = tt.size();

  if(lt < lb) {
    bb.resize(lt, true);
  } else {
    tt.resize(lb, true);
  }

  return make_pair(alessthanb ? tt : bb, alessthanb ? bb : tt);
} 

struct random01 {
    boost::mt19937 &_state;
    double operator()() {
        static boost::uniform_01<boost::mt19937> rng(_state);
        return rng();
    }

    random01(boost::mt19937 &state) : _state(state) {}
};

struct randomMultinomial {
    boost::mt19937 &_state;
        
    template<typename Iter>
    int operator()(Iter first, Iter last) {
        static boost::random::discrete_distribution<> rng(first, last);
        return rng(_state)+1;
    }

    randomMultinomial(boost::mt19937 &state) : _state(state) {}
};

#endif