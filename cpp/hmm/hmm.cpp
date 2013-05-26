#include <iostream>
#include <process.h>

#include <boost/numeric/ublas/matrix_proxy.hpp>

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

/*
void HMM::setKeys(vector<uint> keys) {
  this->keys = keys;
}
*/

// Stablish K keys for N states (equally for each state).
// i.e. K = 9, N = 3 => keys = {1, 1, 1, 2, 2, 2, 3, 3, 3};
vector<uint> HMM::generateKeys(uint N, uint K) {
  vector<uint>keys(K);
  
  uint KN = K / N;

  uint l = KN * N;
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

  return keys;
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

// Simulate T samples for a HMM with parameters p
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

void HMM::EM(params &p, uint MAX_ITER_HMM) {
  if(MAX_ITER_HMM != 0) this->MAX_ITER_HMM = MAX_ITER_HMM;

  for(uint i=0; i<MAX_ITER_HMM; ++i) {

  }
}

bool HMM::ConvergedEM(double ll1, double ll0, double threshold, bool HasIncresed) {
  double delta = fabs(ll1 - ll0);
  double averg  = (fabs(ll1) + fabs(ll0) + DBL_EPSILON) / 2;

  if(HasIncresed) {
    if(delta > 0.5*fabs(ll1)) {
      return true;
    }
  }

  return (delta / averg < threshold);  
}

double HMM::CalculateValues(params p, vector<uint> data, params &q) {
  q = params(p.N, p.K);



  return 0.0;
}

double HMM::BackwardForward(params p, vector<uint> data, 
                            matrix<double> &alpha, matrix <double> &beta,
                            matrix<double> &gamma, matrix <double> &xi) {

  uint N = p.N;
  uint K = p.K;
  uint T = data.size();

  alpha = matrix<double>(N, T);
  beta  = matrix<double>(N, T);
  gamma = matrix<double>(N, T);
  xi    = matrix<double>(N, N);

	double sum, llk;
	double aa, bb;

	vector<double> cn = vector<double>(T);
	matrix<double> xx = matrix<double>(N, N);

	///////////// Forward step /////////////
	for(int n=0; n<N; ++n) {
		alpha(n, 0) = p.priori(0, n) * p.memisn(n, (int)data(0)-1);
	}

  auto ca = column(alpha, 0);
	cn(0) = normalize_row(ca);

	for(int t=1; t<T; ++t) {		
		for(int m=0; m<N; ++m) {
			sum = 0.0;
			for(int n=0; n<N; ++n) {
				aa = alpha(n, t-1) * p.memisn(m, (int)data(t)-1);			
				sum += aa * p.mtrans(n, m);
			}
			alpha(m, t) = sum;
		}

    ca = column(alpha, t);
		cn(t) = normalize_row(ca);
	}

	///////////// Backward step /////////////
	for(int n=0; n<N; ++n) {
		beta(n, T-1) = 1.0;
	}
	
	for(int t=T-2; t>=0; --t) {
		xx *= 0.0;

		for(int m=0; m<N; ++m) {
			sum = 0.0;
			for(int n=0; n<N; ++n) {	
				bb = beta(n, t+1) * p.memisn(n, (int)data(t + 1)-1);			
				sum += bb * p.mtrans(m, n);

				xx(m, n) += p.mtrans(m, n) * alpha(m, t) * bb;
			}
			beta(m, t) = sum;
		}
		//xx.Norm();
    norm2sum(xx);
		xi = xi + xx;

    auto bc = column(beta, t);
		normalize_row(bc);
	}

	// xi.Norm();
  norm2sum(xi);

	///////////// Log Propability /////////////
	sum = 0.0;
	for(int t=0; t<T; ++t) {
		sum += log(cn(t));
	}
	llk = sum;

	///////////// Gamma calc /////////////
	for(int t=0; t<T; ++t) {		
		for(int n=0; n<N; ++n) {
			gamma(n, t) = alpha(n, t) * beta(n, t);
		}

    auto gc = column(gamma, t);
		normalize_row(gc);
	}

	return llk;	
}
