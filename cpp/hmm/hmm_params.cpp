#include "hmm.h"
#include "functions.h"

#include <vector>
#include <set>
#include <algorithm>

#include <process.h>
#include <boost/random.hpp> 

HMM::params::params(uint N, uint K, bool rand) {  
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

matrix<double> HMM::params::estimPriori(vector<uint> data) {
  int N = *(max_element(data.begin(), data.end()));

  matrix<double> priori = matrix<double>(1, N);
  priori(0, data(0) - 1);

  return priori;
}

matrix<double> HMM::params::estimMEmisn(vector<uint> data, vector<uint> keys) {
  auto vec = resize(data, keys);

  data = vec.first;
  keys = vec.second;

  std::set<uint> sata(data.begin(), data.end());
  std::set<uint> seys(keys.begin(), keys.end());

  uint T = data.size();
  uint N = sata.size();
  uint K = seys.size();

  matrix<double> memisn = zero_matrix<double>(N, K);
  
  for(uint i=0; i<N; ++i) {
    auto it = data.begin(), end = data.end();
    int j = 0;
    while(it != end) {
       memisn(*(it++), keys(j++))++;
     }
  }

  normalize(memisn);

  return memisn;
}