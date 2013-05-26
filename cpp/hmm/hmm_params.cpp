#include "hmm.h"
#include "functions.h"

#include <vector>
#include <set>
#include <algorithm>

#include <process.h>
#include <boost/random.hpp> 

HMM::params::params(uint N, uint K, bool rand) : N(N), K(K) {  
  using namespace boost::numeric::ublas;
  using namespace boost::random;
  
  priori = zero_matrix<double>(1, N);
  mtrans = zero_matrix<double>(N, N);
  memisn = zero_matrix<double>(N, K);
    
  if(!rand) return;

  static mt19937 rng(std::time(NULL) + _getpid());
  random01 rnd01(rng);

  std::generate(priori.data().begin(), priori.data().end(), rnd01);
  std::generate(mtrans.data().begin(), mtrans.data().end(), rnd01);
  std::generate(memisn.data().begin(), memisn.data().end(), rnd01);

  normalize(priori);
  normalize(mtrans);
  normalize(memisn);
}

std::ostream& operator<< (std::ostream &out, HMM::params &p) {
  out << p.priori << std:: endl
      << p.mtrans << std:: endl
      << p.memisn << std:: endl;

  return out;
}

matrix<double> HMM::params::estimPriori(vector<uint> observ) {
  int N = *(max_element(observ.begin(), observ.end()));

  matrix<double> priori = matrix<double>(1, N);
  priori(0, observ(0) - 1);

  return priori;
}

matrix<double> HMM::params::estimMEmisn(vector<uint> observ, vector<uint> hidden) {
  auto vec = resize(observ, hidden);

  observ = vec.first;
  hidden = vec.second;

  std::set<uint> sbserv(observ.begin(), observ.end());
  std::set<uint> sidden(hidden.begin(), hidden.end());

  uint T = observ.size();
  uint N = sbserv.size();
  uint K = sidden.size();

  matrix<double> memisn = zero_matrix<double>(N, K);
  
  for(uint i=0; i<N; ++i) {
    auto it = observ.begin(), end = observ.end();
    int j = 0;
    while(it != end) {
       memisn(*(it++), hidden(j++))++;
     }
  }

  normalize(memisn);

  return memisn;
}

vector<uint> estimHidden(vector<uint> observ, vector<uint> keys) {
  int T = observ.size();

  vector<uint> hidden(T);

  for(int t=0; t<T; ++t) {
    hidden(t) = keys(observ(t));
  }

  return hidden;
}