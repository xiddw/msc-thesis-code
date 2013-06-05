#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>

#include <fstream>
#include <utility>

typedef unsigned int uint;
using namespace boost::numeric::ublas;

template class matrix<double>;

// Normalize a boost::vector to sum up to one
/// TODO: check if...
template<typename V>
double normalize_vector(V &v) { 
   double a = sum(v);
   a += (a == 0.0);
   v /= a;
   return a;
} 

// Normalize a boost::matrix in order to have its rows/columns to sum up to one
template<typename M>
vector<double> normalize(M &m, bool byrow = true) { 
  vector<double> s(m.size1());

  for(uint i = 0; i<m.size1(); ++i) { 
     auto r = row(m, i);
     s(i) = normalize_vector(r);
   }

  return s;
} 

// Normalize a boost::matrix in order to have its rows/columns to sum up to one
template<typename M>
double norm2sum(M &m) { 
  //using namespace boost::numeric::ublas;
  
  double s = 0.0;
  
  for(uint i = 0; i<m.size1(); ++i) { 
    for(uint j = 0; j<m.size2(); ++j) { 
     s += m(i, j);
    }
  }

	s += (s == 0.0);

  for(uint i = 0; i<m.size1(); ++i) { 
    for(uint j = 0; j<m.size2(); ++j) { 
     m(i, j) /= s;
    }
  }
  return s;
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

template<typename U>
vector<U> readCSV(std::string file) {
  std::ifstream csvfile;
  csvfile.open(file);

  std::vector<uint> data;

  std::string value;
  while(getline(csvfile, value)) {
    data.push_back(std::stoi(value));
  }     

  csvfile.close();

  uint T = data.size();
  vector<U> res(T);

  for(int t=0; t<T; ++t) {
    res(t) = data.at(t);
  }

  return res;
}

template<typename U>
void writeCSV(vector<U> data, std::string file) {
  std::ofstream csvfile;
  csvfile.open(file);

  uint T = data.size();

  for(int t=0; t<T-1; ++t) {
    csvfile << data(t) << ", ";
  }
  csvfile << data(T-1);
  csvfile.close();
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