#include <iostream>
#include <process.h>

#include "hmm.h"
#include "functions.h"

const uint HMM::DEF_MAX_ITER_ESTIM = 30;
const uint HMM::DEF_MAX_ITER_HMM = 300;

HMM::HMM(uint K, uint T, uint EX, uint MAX_ITER_ESTIM, uint MAX_ITER_HMM) {
  this->K = K;
  this->T = T;
  this->EX = EX;

  this->MAX_ITER_ESTIM = MAX_ITER_ESTIM;
  this->MAX_ITER_HMM = MAX_ITER_HMM;
}

HMM::HMM(vector<uint> data,
         uint MAX_ITER_ESTIM, uint MAX_ITER_HMM) {
  this->data = data;

  K = *(max_element(data.begin(), data.end()));
  T = this->data.size();
  EX = 1;
  
  this->MAX_ITER_ESTIM = MAX_ITER_ESTIM;
  this->MAX_ITER_HMM = MAX_ITER_HMM;
}

void HMM::setKeys(vector<uint> keys) {
  this->keys = keys;
}

// Stablish K keys for N states (equally for each state).
// i.e. K = 9, N = 3 => keys = {1, 1, 1, 2, 2, 2, 3, 3, 3};
void HMM::genKeys(uint N) {
  this->N = N;
  keys = vector<uint>(K);
  
  uint KN = K / this->N;

  uint l = KN * this->N;
  int k = 1;

  for(uint i=0; i<l; ++i) {    
    if(i % KN == 0 && i>0) k++;
    keys[i] = k;
  }

  if(l < K) {
    for(uint i=l; i<K; i++) {
      keys[i] = k;
    }
  }
}

// Generate parameters for an HMM model 
// a priori matrix, (1, 0, ..., 0)
// a transition matrix (diagonally dominant) and
// a emision matrix (diagonally dominant by blocks)
HMM::params HMM::generateParams(uint N, uint K) {
  HMM::params p(N, K);

  p.priori(0, 0) = 1;

  uint a = 1;
  uint b = 50;
  uint c = 5;

  p.mtrans(0, 0) = b;
  for(uint i=1; i<N; ++i) {
    p.mtrans(i, i) = b;
    p.mtrans(i-1, i) = a;
    p.mtrans(i, i-1) = c;
  }

  uint KN = K / N;

  uint l = KN * N;
  int k = 0;

  for(uint i=0; i<l; ++i) {    
    if(i % KN == 0 && i>0) k++;
    p.memisn(k, i) = 100;
  }

  if(l < K) {
    for(uint i=l; i<K; i++) {
      p.memisn(k, i) = 100;
    }
  }

  normalize(p.priori);
  normalize(p.mtrans);
  normalize(p.memisn);

  return p;
}

vector<uint> HMM::sample(params p, uint T) {
  vector<uint> data(T);

  static boost::mt19937 rng(std::time(NULL) + _getpid());
  randomMultinomial rndMult(rng);

  uint spkr, word;
  auto row_priori = row(p.priori, 0);

  spkr = rndMult(row_priori.begin(), row_priori.end());

  for(int t=0; t<T; ++t) {
    auto row_memisn = row(p.memisn, spkr-1);
    auto row_mtrans = row(p.mtrans, spkr-1);

    word = rndMult(row_memisn.begin(), row_memisn.end());
    spkr = rndMult(row_mtrans.begin(), row_mtrans.end());

    data(t) = spkr;
  }

  return data;
}