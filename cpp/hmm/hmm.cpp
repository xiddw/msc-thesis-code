#include <iostream>
#include <process.h>

#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "hmm.h"
#include "functions.h"

const uint HMM::DEF_MAX_ITER_ESTIM = 350;
const uint HMM::DEF_MAX_ITER_HMM = 30;

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

  for(uint t=0; t<T; ++t) {
    auto row_memisn = row(p.memisn, spkr-1);
    auto row_mtrans = row(p.mtrans, spkr-1);

    word = rndMult(row_memisn.begin(), row_memisn.end());
    spkr = rndMult(row_mtrans.begin(), row_mtrans.end());

    data(t) = spkr;
  }

  return data;
}

double HMM::EM(params &p, uint MAX_ITER_ESTIM) {
  // if MAX_ITER_HMM != 0, then assing it. Otherwise, use constructor (default) value
  if(MAX_ITER_ESTIM != 0) this->MAX_ITER_ESTIM = MAX_ITER_ESTIM;

  double ll0, ll1;
  ll0 = -LONG_MAX;
  ll1 = 0.0;

  // vector<double> LL = vector<double>(this->MAX_ITER_ESTIM);
  matrix<double> gamma;

  std::cout << "total iter: " << this->MAX_ITER_ESTIM << std::endl;

  for(uint i=0; i< this->MAX_ITER_ESTIM; ++i) {
    params q = p;

    ll1 = CalculateValues(p, data, q, gamma);

    //if(ConvergedEM(ll0, ll1)) break;

    normalize(q.priori);
    normalize(q.mtrans);
    normalize(q.memisn);

    p = q;
    // ll0 = ll1;
    // LL(i) = ll0;

    // std::cout << "iter: " << i << std::endl;
    // std::cout << "ll(" << i << "): " << LL(i) << std::endl;
    std::cout << ".";
  }
  std::cout << std::endl << "ll: " << ll1 << std::endl;
  // std::cout << std::endl << "ll: " << LL(this->MAX_ITER_ESTIM-1) << std::endl;

  p.hidden = vector<uint>(T);
  for(uint t=0; t<T; ++t) {
    auto cg = column(gamma, t);

    uint k = distance(cg.begin(), max_element(cg.begin(), cg.end()));
    p.hidden(t) = k;
  }

  return ll1;
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

double HMM::CalculateValues(params p, vector<uint> data, params &q, matrix<double> &gamma) {
  q = params(p.N, p.K);

  uint N = p.N;
  uint K = p.K;
  uint T = data.size();

  double ll = 0.0;

  // if it would be multiple series then cycle over them
  matrix<double> alpha, beta, xi;
  
  ll = BackwardForward(p, data, alpha, beta, gamma, xi);

  // update priori matrix
  auto g1 = column(gamma, 0);
  std::copy(g1.begin(), g1.end(), q.priori.begin2());

  // update transition matrix
  q.mtrans = xi;

  // update emision matrix
  uint k;
  for(uint t=0; t<T; ++t) {
    k = data(t) - 1;
    for(uint n=0; n<N; ++n) {
      q.memisn(n, k) += gamma(n, t);
    }
  }

	// std::cout << "logl:" << std::endl;
	// std::cout << ll << std::endl;

  return ll;
}

double HMM::BackwardForward(params p, vector<uint> data, 
                            matrix<double> &alpha, matrix <double> &beta,
                            matrix<double> &gamma, matrix <double> &xi) {

  uint N = p.N;
  uint K = p.K;
  uint T = data.size();

  alpha = zero_matrix<double>(N, T);
  beta  = zero_matrix<double>(N, T);
  gamma = zero_matrix<double>(N, T);
  xi    = zero_matrix<double>(N, N);

	double sum, llk;
	double aa, bb;

	vector<double> cn = zero_vector<double>(T);
	matrix<double> xx;

	///////////// Forward step /////////////
	for(uint n=0; n<N; ++n) {
		alpha(n, 0) = p.priori(0, n) * p.memisn(n, (int)data(0)-1);
	}
	
	auto ac = column(alpha, 0);
	cn(0) = normalize_vector(ac);

	for(uint t=1; t<T; ++t) {		
		for(uint m=0; m<N; ++m) {
			sum = 0.0;
			for(uint n=0; n<N; ++n) {
				aa = alpha(n, t-1) * p.memisn(m, (int)data(t)-1);			
				sum += aa * p.mtrans(n, m);
			}
			alpha(m, t) = sum;
		}

    auto ac = column(alpha, t);
		cn(t) = normalize_vector(ac);
	}

	///////////// Backward step /////////////
	for(uint n=0; n<N; ++n) {
		beta(n, T-1) = 1.0;
	}
	
	for(int t=T-2; t>=0; --t) {
		xx = zero_matrix<double>(N, N);

		for(uint m=0; m<N; ++m) {
			sum = 0.0;
			for(uint n=0; n<N; ++n) {	
				bb = beta(n, t+1) * p.memisn(n, (int)data(t+1)-1);			
				sum += bb * p.mtrans(m, n);

				xx(m, n) += p.mtrans(m, n) * alpha(m, t) * bb;
			}
			beta(m, t) = sum;
		}
    norm2sum(xx); 
		xi = xi + xx;

    auto bc = column(beta, t);
		normalize_vector(bc);
	}
  
  norm2sum(xi);

  // std::cout << "CN:" << std::endl;
  // std::cout << cn << std::endl;

	///////////// Log Propability /////////////
	sum = 0.0;
	for(uint t=0; t<T; ++t) {
		sum += log(cn(t));
	}
	llk = sum;

	///////////// Gamma calc /////////////
	for(uint t=0; t<T; ++t) {		
		for(uint n=0; n<N; ++n) {
			gamma(n, t) = alpha(n, t) * beta(n, t);
		}

    auto gc = column(gamma, t);
		normalize_vector(gc);
	}

	return llk;	
}
