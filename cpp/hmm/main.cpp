#include <iostream>
#include <process.h>

#include <process.h>
#include <boost/random.hpp> 

#include "hmm.h"
#include "functions.h"

int main() {

  auto data = readCSV<uint>("E:\\ESCUELA\\CIMAT\\4 Semestre\\ST2\\prog\\pruebas\\cuervo1f_120.csv");

  uint nn = 6;
  uint kk = 120;    
  uint tt = data.size();

  HMM::params p(nn, kk, true);

  HMM h(data);
  h.EM(p);

  auto res = p.hidden;
  writeCSV<uint>(res, "E:\\ESCUELA\\CIMAT\\4 Semestre\\ST2\\prog\\pruebas\\cuervo1f_mground.csv");

  getchar();
  return 0;
}