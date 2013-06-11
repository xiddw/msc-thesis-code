#include <iostream>
#include <process.h>

#include <process.h>
#include <boost/random.hpp> 

#include "hmm.h"
#include "functions.h"

int main() {
  /*
	using namespace boost::numeric::ublas;
  matrix<double> m (3, 4);
  for (unsigned i = 0; i < m.size1 (); ++ i)
      for (unsigned j = 0; j < m.size2 (); ++ j)
          m (i, j) = 4 * i + j;

  std::cout << "original: " << std::endl;
  std::cout << m << std::endl;

  std::cout << "range: " << std::endl;
  matrix_range<matrix<double> > g1(m, range(0, 0), range(0, 4));
  std::cout << g1 << std::endl;

  std::cout << "row: " << std::endl;
  auto mr = row(m, 1);
  matrix<double> g2(4, 1);
  copy(mr.begin(), mr.end(), g2.begin1());
  std::cout << g2 << std::endl;

  getchar();
  
  auto rr = column(m, 2);
  std::cout << "RR: " << std::endl;
  std::cout << rr << std::endl;
  std::cout << "sum(RR): " << std::endl;
  std::cout << sum(rr) << std::endl;
  std::cout << "norm(RR): " << std::endl;
  std::cout << rr / sum(rr) << std::endl;
  std::cout << "norm(RR): " << std::endl;
  std::cout << normalize_vector(rr) << std::endl;
  std::cout << "RR2: " << std::endl;
  std::cout << rr << std::endl;
  //*/

  //std::string path = "C:\\Users\\Estudiante\\Documents\\GitHub\\msc-thesis-code\\";
  std::string path = "E:\\ESCUELA\\CIMAT\\4 Semestre\\ST2\\prog\\";
	
  auto data = readCSV<uint>(path + "pruebas\\cuervo1f_120.csv");

	//data.resize(4, true);

  uint nn = 6;
  uint kk = 120;    
  uint tt = data.size();

  HMM::params m(nn, kk);

  HMM h(data);
  
  double ll = -LONG_MAX, mm = -LONG_MAX;
  uint ii = 0;
  for(uint i=0; i<30; ++i) {
    HMM::params p(nn, kk, true);
    ll = h.EM(p);

    if(ll > mm) {
      ii = i;
      mm = ll;
      m = p;
    }
  }  

  std::cout << "best overal (iter=>" << ii << "): " << mm << std::endl;

  auto res = m.hidden;
  writeCSV<uint>(res, path + "pruebas\\cuervo1f_mground.csv");

  getchar();
  return 0;
}