#include <iostream>
#include <process.h>

#include <process.h>
#include <boost/random.hpp> 

#include <ctime>
#include <omp.h>

#include "hmm.h"
#include "functions.h"



int main() {
  std::string path = "C:\\Users\\Estudiante\\Documents\\GitHub\\msc-thesis-code\\";
  // std::string path = "E:\\ESCUELA\\CIMAT\\4 Semestre\\ST2\\prog\\";
	
  auto data = readCSV<uint>(path + "pruebas\\cuervo1f_120.csv");

  uint nn = 6;
  uint kk = 120;    
  uint tt = data.size();

  HMM::params m(nn, kk);

  HMM h(data);
  
	int max = h.MaxIterHMM();

	int i, tid, nthreads;

	// Auxiliar variables for omp process :)
	double *ll1 = new double[max], 
				 *ll2 = new double[max], 
					mm1, mm2;

	HMM::params *pp1 = new HMM::params[max], 
							*pp2 = new HMM::params[max], 
							 oo1, oo2;

	clock_t ts = clock();

	#pragma omp parallel for private(tid, i) shared(ll1, ll2, pp1, pp2) 
	for(i=0; i<max; ++i) {
		HMM::params p(nn, kk, true);
		ll1[i] = h.EM(p);
		pp1[i] = p;

		HMM::params q(nn+1, kk, true);
		ll2[i] = h.EM(q);
		pp2[i] = q;

		std::cout << "."; 
	}

	mm1 = ll1[0];
	mm2 = ll2[0];

	int ii1 = 0, ii2;
	for(i=1; i<max; ++i) {
		if(ll1[i] > mm1) {
			mm1 = ll1[i];
			ii1 = i;
		}

		if(ll2[i] > mm2) {
			mm2 = ll2[i];
			ii2 = i;
		}
	}
	oo1 = pp1[ii1];
	oo2 = pp2[ii2];

	clock_t te = clock();

	printf("\n Time taken: %.2fs\n", (double)(te -ts)/CLOCKS_PER_SEC);

  std::cout << "best overal [model 1] (iter=>" << ii1 << "): " << mm1 << std::endl;
	std::cout << "best overal [model 2] (iter=>" << ii2 << "): " << mm2 << std::endl;

  auto res1 = oo1.hidden;
	auto res2 = oo2.hidden;
  writeCSV<uint>(res1, path + "pruebas\\cuervo1f_" + std::to_string(kk) + "mground.csv");
	writeCSV<uint>(res2, path + "pruebas\\cuervo1f_" + std::to_string(kk+1) + "mground.csv");

  getchar();
  return 0;
}